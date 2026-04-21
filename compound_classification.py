"""
compound_classification.py - Step 2c of the GCMS processing pipeline.

Queries PubChem for each identified compound (non-empty compound_name in
feature_metadata.csv) to retrieve the canonical SMILES and InChIKey, then
classifies the structure using two independent web services:

  ClassyFire (classyfire.wishartlab.com)
    - Kingdom, Superclass, Class, Subclass, Direct Parent
    - Rule-based chemical taxonomy; good coverage for terpenoids, lipids,
      aldehydes, alcohols, and other GC-MS volatiles
    - Queried by InChIKey via the entities endpoint (direct API, higher
      coverage than the PubChem PUG View mirror)

  NPClassifier (npclassifier.gnps2.org)
    - Pathway, Superclass, Class
    - Deep-learning natural-product classifier; provides biogenic pathway
      information (e.g. Terpenoids → Sesquiterpenoids → Bisabolane sesqui-
      terpenoids).  Best suited for plant/microbial secondary metabolites.
    - Queried by canonical SMILES

API flow per compound:
  1. PubChem name  → CID                      (PUG REST)
  2. PubChem CID   → IUPACName + MolecularFormula  (PUG REST)
  3. PubChem CID   → canonical SMILES + InChIKey   (PUG REST)
  4. ClassyFire InChIKey → taxonomy           (classyfire.wishartlab.com)
  5. NPClassifier SMILES → NP class           (npclassifier.gnps2.org)

All results are cached in output/pubchem_cache.json so re-runs skip
already-fetched data without making additional network requests.
HTTP 503 responses are retried with exponential backoff on all endpoints.

Rate limits honoured:
  PubChem      < 5 req/s  (PUBCHEM_RATE_LIMIT_DELAY,   default 0.35 s)
  ClassyFire   ~1 req/s   (CLASSYFIRE_RATE_LIMIT_DELAY, default 1.0 s;
                            no stated limit — conservative by default)
  NPClassifier ~1 req/s   (NPCLASSIFIER_RATE_LIMIT_DELAY, default 1.0 s;
                            no stated limit — conservative by default)

API usage policies:
  PubChem     : https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access
  ClassyFire  : https://classyfire.wishartlab.com  (academic use)
  NPClassifier: https://github.com/mwang87/NP-Classifier

Input  : output/feature_metadata.csv
Output : output/compound_classes.csv
         output/feature_metadata_enriched.csv
         output/pubchem_cache.json  (cache; persists across runs)

Usage:
    python compound_classification.py
"""

import json
import os
import sys
import time

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import pandas as pd
import requests

import config


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_BASE_PUG          = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_BASE_CLASSYFIRE   = "http://classyfire.wishartlab.com"
_BASE_NPCLASSIFIER = "https://npclassifier.gnps2.org"

# Sentinel stored in the cache to mark a name / CID as permanently not found.
# This prevents re-querying entries that are known to have no record.
_NOT_FOUND = "__NOT_FOUND__"

# Bump this string whenever the ClassyFire or NPClassifier fetch logic changes
# in a way that may produce different results from cached data.
# On first run after a version bump, stale _NOT_FOUND entries in
# cid_to_classyfire and cid_to_npclassifier are cleared so they can be
# re-tried with the updated method.
_CLASSYFIRE_API_VERSION = "direct_v1"

# Classification status labels written to the output CSV.
STATUS_FOUND       = "found"
STATUS_CID_MISSING = "cid_not_found"
STATUS_CF_MISSING  = "classyfire_not_available"
STATUS_UNNAMED     = "unnamed"
STATUS_ERROR       = "error"


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def _load_cache(cache_path):
    """Load the JSON cache file; return a safe defaulted dict on any error."""
    empty = {
        "name_to_cid":         {},
        "cid_to_properties":   {},
        "cid_to_smiles":       {},
        "cid_to_classyfire":   {},
        "cid_to_npclassifier": {},
    }
    if not os.path.exists(cache_path):
        return empty
    try:
        with open(cache_path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        # Ensure all expected keys exist (handles caches from older versions)
        for key in empty:
            data.setdefault(key, {})
        return data
    except (json.JSONDecodeError, OSError) as exc:
        print(f"  [classification] WARNING: could not read cache ({exc}); starting fresh.")
        return empty


def _save_cache(cache, cache_path):
    """Write-through: persist the cache immediately after every update."""
    try:
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, indent=2)
    except OSError as exc:
        print(f"  [classification] WARNING: could not write cache ({exc}).")


# ---------------------------------------------------------------------------
# HTTP helper
# ---------------------------------------------------------------------------

def _get_with_retry(url, headers, rate_delay, max_retries=5):
    """
    Sleep rate_delay seconds, then GET the URL.
    Retries on HTTP 503 with exponential back-off.

    Returns
    -------
    (response, is_not_found) where:
      - response is None on unrecoverable error
      - is_not_found is True for HTTP 400 or 404
        (400: no data for this heading / bad request;
         404: record does not exist)
    """
    time.sleep(rate_delay)          # always honour the rate limit first
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, headers=headers, timeout=20)
        except requests.RequestException as exc:
            wait = 2 ** attempt
            print(f"  [classification] network error ({exc}); retrying in {wait}s …")
            time.sleep(wait)
            continue

        if resp.status_code == 200:
            return resp, False
        if resp.status_code in (400, 404):
            # 400: no data for this heading (PubChem / ClassyFire) or bad request
            # 404: record not found
            return None, True
        if resp.status_code in (429, 503):
            wait = 2 ** attempt
            print(f"  [classification] {resp.status_code} received; backing off {wait}s …")
            time.sleep(wait)
            continue
        # Any other unexpected status: log and abort retries
        print(f"  [classification] unexpected HTTP {resp.status_code} for {url}")
        return None, False

    print(f"  [classification] gave up after {max_retries} attempts: {url}")
    return None, False


# ---------------------------------------------------------------------------
# PubChem API fetch functions
# ---------------------------------------------------------------------------

def _fetch_cid(name, headers, rate_delay):
    """
    Look up a compound name in PubChem.
    Returns the first CID as a string, or None if not found / error.
    """
    encoded = requests.utils.quote(name, safe="")
    url     = f"{_BASE_PUG}/compound/name/{encoded}/cids/JSON"
    resp, is_404 = _get_with_retry(url, headers, rate_delay)
    if is_404 or resp is None:
        return None
    try:
        cids = resp.json()["IdentifierList"]["CID"]
        return str(cids[0])
    except (KeyError, IndexError, ValueError):
        return None


def _fetch_properties(cid, headers, rate_delay):
    """
    Fetch MolecularFormula and IUPACName for a CID.
    Returns dict with keys molecular_formula, iupac_name, or None on failure.
    """
    url  = f"{_BASE_PUG}/compound/cid/{cid}/property/MolecularFormula,IUPACName/JSON"
    resp, is_404 = _get_with_retry(url, headers, rate_delay)
    if is_404 or resp is None:
        return None
    try:
        props = resp.json()["PropertyTable"]["Properties"][0]
        return {
            "molecular_formula": props.get("MolecularFormula"),
            "iupac_name":        props.get("IUPACName"),
        }
    except (KeyError, IndexError, ValueError):
        return None


def _fetch_smiles_inchikey(cid, headers, rate_delay):
    """
    Fetch canonical SMILES and InChIKey for a CID from PubChem.
    Returns dict with keys smiles, inchikey, or None on failure.

    PubChem returns different SMILES keys depending on the compound:
      CanonicalSMILES    - most compounds with defined stereochemistry
      IsomericSMILES     - stereochemistry-aware variant
      ConnectivitySMILES - compounds with undefined stereocentres (e.g. many
                           GC-MS library matches); CID 11644 (chloroiodomethane)
                           is a known example where only this key is returned
    All three are requested; the first non-empty value is used.
    """
    url  = (f"{_BASE_PUG}/compound/cid/{cid}/property/"
            f"CanonicalSMILES,IsomericSMILES,ConnectivitySMILES,InChIKey/JSON")
    resp, is_404 = _get_with_retry(url, headers, rate_delay)
    if is_404 or resp is None:
        return None
    try:
        props = resp.json()["PropertyTable"]["Properties"][0]
        smiles = (props.get("CanonicalSMILES")
                  or props.get("IsomericSMILES")
                  or props.get("ConnectivitySMILES"))
        return {
            "smiles":   smiles,
            "inchikey": props.get("InChIKey"),
        }
    except (KeyError, IndexError, ValueError):
        return None


# ---------------------------------------------------------------------------
# ClassyFire direct API
# ---------------------------------------------------------------------------

def _fetch_classyfire_direct(inchikey, headers, rate_delay):
    """
    Look up a compound in ClassyFire by InChIKey.

    Uses the ClassyFire entities endpoint (classyfire.wishartlab.com).
    This direct API has substantially higher coverage than the PubChem
    PUG View ClassyFire mirror, which is incomplete for many GC-MS volatiles.

    Returns a dict with keys kingdom, superclass, class_, subclass,
    direct_parent (all may be None), or None if not found / no data.

    API policy: academic / non-commercial use.
    https://classyfire.wishartlab.com
    """
    url  = f"{_BASE_CLASSYFIRE}/entities/{inchikey}.json"
    resp, is_404 = _get_with_retry(url, headers, rate_delay)
    if is_404 or resp is None:
        return None
    try:
        data = resp.json()
    except ValueError:
        return None

    def _name(obj):
        return obj.get("name") if isinstance(obj, dict) else None

    kingdom       = _name(data.get("kingdom"))
    superclass    = _name(data.get("superclass"))
    class_        = _name(data.get("class"))
    subclass      = _name(data.get("subclass"))
    direct_parent = _name(data.get("direct_parent"))

    if not any([kingdom, superclass, class_, subclass, direct_parent]):
        return None

    return {
        "kingdom":       kingdom,
        "superclass":    superclass,
        "class_":        class_,
        "subclass":      subclass,
        "direct_parent": direct_parent,
    }


# ---------------------------------------------------------------------------
# NPClassifier direct API
# ---------------------------------------------------------------------------

def _fetch_npclassifier_direct(smiles, headers, rate_delay):
    """
    Classify a compound by canonical SMILES using the NPClassifier API.

    Uses npclassifier.gnps2.org. Returns a dict with keys
    npclassifier_pathway, npclassifier_superclass, npclassifier_class,
    or None if no classification is available.

    NPClassifier is optimised for plant and microbial secondary metabolites
    (terpenoids, alkaloids, polyketides, shikimates). Coverage for simple
    industrial volatiles and organohalogens is lower.

    API policy: public tool (UCSD/GNPS project); no stated rate limit.
    https://github.com/mwang87/NP-Classifier
    """
    encoded = requests.utils.quote(smiles, safe="")
    url     = f"{_BASE_NPCLASSIFIER}/classify?smiles={encoded}"
    resp, is_404 = _get_with_retry(url, headers, rate_delay)
    if is_404 or resp is None:
        return None
    try:
        data = resp.json()
    except ValueError:
        return None

    pathway    = data.get("pathway_results",    [])
    superclass = data.get("superclass_results", [])
    class_     = data.get("class_results",      [])

    if not any([pathway, superclass, class_]):
        return None

    return {
        "npclassifier_pathway":    pathway[0]    if pathway    else None,
        "npclassifier_superclass": superclass[0] if superclass else None,
        "npclassifier_class":      class_[0]     if class_     else None,
    }


# ---------------------------------------------------------------------------
# Name filter
# ---------------------------------------------------------------------------

def _is_unnamed(name):
    """
    Return True for names that should not be queried:
      - NaN / None / non-string
      - empty or whitespace-only
      - TraceFinder auto-generated placeholders such as "Peak@3.102"
    """
    if not isinstance(name, str):
        return True
    name = name.strip()
    if not name:
        return True
    if name.lower().startswith("peak@"):
        return True
    return False


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run(cfg=config):
    """Query PubChem + ClassyFire + NPClassifier for compound annotations."""

    # --- Config parameters --------------------------------------------------
    user_agent    = getattr(cfg, "PUBCHEM_USER_AGENT",
                            "TF_NT_pipeline/1.0 (nontargeted GCMS metabolomics; "
                            "contact: user@example.com)")
    pubchem_delay = getattr(cfg, "PUBCHEM_RATE_LIMIT_DELAY", 0.35)
    cf_delay      = getattr(cfg, "CLASSYFIRE_RATE_LIMIT_DELAY", 1.0)
    npc_delay     = getattr(cfg, "NPCLASSIFIER_RATE_LIMIT_DELAY", 1.0)
    cache_path    = getattr(cfg, "PUBCHEM_CACHE_FILE",
                            os.path.join(cfg.OUTPUT_DIR, "pubchem_cache.json"))
    cache_only    = getattr(cfg, "PUBCHEM_CACHE_ONLY", False)

    meta_path    = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")
    out_classes  = os.path.join(cfg.OUTPUT_DIR, "compound_classes.csv")
    out_enriched = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")

    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 2c: compound classification --------------------------------")
    print(f"   feature_metadata  : {meta_path}")
    print(f"   cache file        : {cache_path}")
    if cache_only:
        print("   mode              : cache-only (no network requests)")
    else:
        print(f"   PubChem delay     : {pubchem_delay} s/request")
        print(f"   ClassyFire delay  : {cf_delay} s/request")
        print(f"   NPClassifier delay: {npc_delay} s/request")

    # --- Load feature metadata ----------------------------------------------
    if not os.path.exists(meta_path):
        sys.exit(
            f"[classification] Input not found: {meta_path}\n"
            "Run data_import.py first."
        )

    meta = pd.read_csv(meta_path, index_col="feature_id")

    if "compound_name" not in meta.columns:
        print("  [classification] No 'compound_name' column found; nothing to classify.")
        return

    cache   = _load_cache(cache_path)
    headers = {"User-Agent": user_agent}

    # --- Cache migration: clear null-SMILES entries -------------------------
    # A previous run stored {"smiles": null, "inchikey": "..."} for all CIDs
    # because the CanonicalSMILES key was silently missing from the response.
    # Clear those so the fetch is retried with the corrected parser.
    null_smiles = [k for k, v in cache["cid_to_smiles"].items()
                   if isinstance(v, dict) and v.get("smiles") is None]
    if null_smiles:
        print(f"   [migration] clearing {len(null_smiles)} null-SMILES cache entries "
              f"— will re-fetch CanonicalSMILES")
        for k in null_smiles:
            del cache["cid_to_smiles"][k]
        # Also clear NPClassifier since it depends on SMILES
        cache["cid_to_npclassifier"] = {}
        _save_cache(cache, cache_path)

    # --- Cache migration: clear stale _NOT_FOUND entries from old API -------
    # Previous versions fetched ClassyFire via PubChem PUG View, which has
    # incomplete coverage and returns HTTP 400 for missing sections.  Those
    # _NOT_FOUND sentinels are cleared here so the direct API gets a chance.
    if cache.get("__classyfire_api_version__") != _CLASSYFIRE_API_VERSION:
        n_cf  = sum(1 for v in cache["cid_to_classyfire"].values()   if v == _NOT_FOUND)
        n_npc = sum(1 for v in cache["cid_to_npclassifier"].values() if v == _NOT_FOUND)
        if n_cf or n_npc:
            print(f"   [migration] clearing {n_cf} ClassyFire + {n_npc} NPClassifier "
                  f"stale cache entries — will re-fetch via direct APIs")
        cache["cid_to_classyfire"]   = {k: v for k, v in cache["cid_to_classyfire"].items()
                                         if v != _NOT_FOUND}
        cache["cid_to_npclassifier"] = {k: v for k, v in cache["cid_to_npclassifier"].items()
                                         if v != _NOT_FOUND}
        cache["__classyfire_api_version__"] = _CLASSYFIRE_API_VERSION
        _save_cache(cache, cache_path)

    named_mask = ~meta["compound_name"].apply(_is_unnamed)
    n_named    = named_mask.sum()
    n_total    = len(meta)

    print(f"   features total    : {n_total}")
    print(f"   named features    : {n_named}  (unnamed / Peak@ are skipped)")

    names_to_query  = meta.loc[named_mask, "compound_name"].unique().tolist()
    new_cid_lookups = sum(1 for n in names_to_query if n not in cache["name_to_cid"])
    if cache_only:
        hits = len(names_to_query) - new_cid_lookups
        print(f"   cache hits        : {hits} / {len(names_to_query)} named compounds")
        print(f"   skipped (no cache): {new_cid_lookups}  (set PUBCHEM_CACHE_ONLY=False to fetch)")
    else:
        print(f"   new CID lookups   : {new_cid_lookups}  (rest served from cache)")

    # --- API query loop -----------------------------------------------------
    for idx, (fid, row) in enumerate(meta[named_mask].iterrows(), start=1):
        name   = row["compound_name"]
        prefix = f"   [{idx}/{n_named}] {name[:50]:<50}"

        # 1. name → CID  (PubChem PUG REST)
        if name not in cache["name_to_cid"]:
            if cache_only:
                continue
            cid = _fetch_cid(name, headers, pubchem_delay)
            cache["name_to_cid"][name] = cid if cid is not None else _NOT_FOUND
            _save_cache(cache, cache_path)

        cid_val = cache["name_to_cid"][name]
        if cid_val == _NOT_FOUND:
            print(f"{prefix}  -> CID not found")
            continue

        cid = cid_val

        # 2. CID → IUPACName + MolecularFormula  (PubChem PUG REST)
        if cid not in cache["cid_to_properties"]:
            if cache_only:
                continue
            props = _fetch_properties(cid, headers, pubchem_delay)
            cache["cid_to_properties"][cid] = props if props is not None else _NOT_FOUND
            _save_cache(cache, cache_path)

        # 3. CID → canonical SMILES + InChIKey  (PubChem PUG REST)
        if cid not in cache["cid_to_smiles"]:
            if cache_only:
                continue
            si = _fetch_smiles_inchikey(cid, headers, pubchem_delay)
            cache["cid_to_smiles"][cid] = si if si is not None else _NOT_FOUND
            _save_cache(cache, cache_path)

        si_val   = cache["cid_to_smiles"].get(cid)
        smiles   = si_val.get("smiles")   if isinstance(si_val, dict) else None
        inchikey = si_val.get("inchikey") if isinstance(si_val, dict) else None

        # 4. InChIKey → ClassyFire taxonomy  (classyfire.wishartlab.com)
        if inchikey and cid not in cache["cid_to_classyfire"]:
            if not cache_only:
                cf = _fetch_classyfire_direct(inchikey, headers, cf_delay)
                cache["cid_to_classyfire"][cid] = cf if cf is not None else _NOT_FOUND
                _save_cache(cache, cache_path)

        # 5. SMILES → NPClassifier  (npclassifier.gnps2.org)
        if smiles and cid not in cache["cid_to_npclassifier"]:
            if not cache_only:
                npc = _fetch_npclassifier_direct(smiles, headers, npc_delay)
                cache["cid_to_npclassifier"][cid] = npc if npc is not None else _NOT_FOUND
                _save_cache(cache, cache_path)

        cf_val  = cache["cid_to_classyfire"].get(cid)
        kingdom = cf_val.get("kingdom") if isinstance(cf_val, dict) else None
        npc_val = cache["cid_to_npclassifier"].get(cid)
        npc_cls = npc_val.get("npclassifier_class") if isinstance(npc_val, dict) else None
        print(f"{prefix}  CID={cid}  kingdom={kingdom or '—'}  NPC={npc_cls or '—'}")

    # --- Build output rows --------------------------------------------------
    rows = []
    for fid, row in meta.iterrows():
        name = row.get("compound_name", "")
        rec  = {
            "feature_id":              fid,
            "compound_name":           name,
            "classification_status":   STATUS_UNNAMED,
            "pubchem_cid":             None,
            "iupac_name":              None,
            "molecular_formula":       None,
            "kingdom":                 None,
            "superclass":              None,
            "class":                   None,
            "subclass":                None,
            "direct_parent":           None,
            "npclassifier_pathway":    None,
            "npclassifier_superclass": None,
            "npclassifier_class":      None,
        }

        if _is_unnamed(name):
            rows.append(rec)
            continue

        cid_val = cache["name_to_cid"].get(name)
        if cid_val is None or cid_val == _NOT_FOUND:
            rec["classification_status"] = STATUS_CID_MISSING
            rows.append(rec)
            continue

        rec["pubchem_cid"] = cid_val

        props = cache["cid_to_properties"].get(cid_val)
        if isinstance(props, dict):
            rec["iupac_name"]        = props.get("iupac_name")
            rec["molecular_formula"] = props.get("molecular_formula")

        cf = cache["cid_to_classyfire"].get(cid_val)
        if isinstance(cf, dict):
            rec["kingdom"]       = cf.get("kingdom")
            rec["superclass"]    = cf.get("superclass")
            rec["class"]         = cf.get("class_")
            rec["subclass"]      = cf.get("subclass")
            rec["direct_parent"] = cf.get("direct_parent")
            rec["classification_status"] = STATUS_FOUND
        else:
            rec["classification_status"] = STATUS_CF_MISSING

        npc = cache["cid_to_npclassifier"].get(cid_val)
        if isinstance(npc, dict):
            rec["npclassifier_pathway"]    = npc.get("npclassifier_pathway")
            rec["npclassifier_superclass"] = npc.get("npclassifier_superclass")
            rec["npclassifier_class"]      = npc.get("npclassifier_class")

        rows.append(rec)

    # --- Write compound_classes.csv -----------------------------------------
    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_classes, index=False)
    print(f"\n   Written: {out_classes}")

    # --- Write feature_metadata_enriched.csv --------------------------------
    merge_cols = [c for c in out_df.columns
                  if c not in ("feature_id", "compound_name")]
    enriched = meta.reset_index().merge(
        out_df[["feature_id"] + merge_cols],
        on="feature_id",
        how="left",
    ).set_index("feature_id")
    enriched.to_csv(out_enriched)
    print(f"   Written: {out_enriched}")

    # --- Write feature_metadata_enriched_targeted.csv ------------------------
    # Filter to only targeted compounds and add delta_rt and delta_mz columns
    targeted_list = getattr(cfg, "TARGETED_LIST", [])
    if not targeted_list:
        targeted_list = getattr(cfg, "EXCLUSION_LIST", [])
    
    if targeted_list:
        # Match targeted compounds by RT and m/z coordinates instead of name
        rt_margin = getattr(cfg, "EXCLUSION_RT_MARGIN", cfg.RT_MARGIN)
        mz_tolerance = getattr(cfg, "EXCLUSION_MZ_TOLERANCE", cfg.MZ_TOLERANCE)
        
        targeted_indices = []
        for entry in targeted_list:
            if len(entry) < 2:
                continue
            target_rt, target_mz = entry[0], entry[1]
            target_name = entry[2] if len(entry) >= 3 else None
            
            # Match features within RT and m/z tolerance
            if target_mz is not None:
                mask = (
                    (enriched["mean_rt"] >= target_rt - rt_margin) &
                    (enriched["mean_rt"] <= target_rt + rt_margin) &
                    (enriched["mean_mz"] >= target_mz - mz_tolerance) &
                    (enriched["mean_mz"] <= target_mz + mz_tolerance)
                )
            else:
                # Match RT only if m/z is None
                mask = (
                    (enriched["mean_rt"] >= target_rt - rt_margin) &
                    (enriched["mean_rt"] <= target_rt + rt_margin)
                )
            
            matches = enriched[mask].index.tolist()
            targeted_indices.extend(matches)
        
        # Remove duplicates while preserving order
        targeted_indices = list(dict.fromkeys(targeted_indices))
        targeted_enriched = enriched.loc[targeted_indices].copy()
        
        # Add delta_rt and delta_mz columns
        targeted_enriched["delta_rt"] = targeted_enriched["rt_max"] - targeted_enriched["rt_min"]
        targeted_enriched["delta_mz"] = targeted_enriched["mz_max"] - targeted_enriched["mz_min"]
        
        # Write to CSV
        out_targeted = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched_targeted.csv")
        targeted_enriched.to_csv(out_targeted)
        print(f"   Written: {out_targeted}")

    # --- Summary ------------------------------------------------------------
    found     = (out_df["classification_status"] == STATUS_FOUND).sum()
    no_cid    = (out_df["classification_status"] == STATUS_CID_MISSING).sum()
    no_cf     = (out_df["classification_status"] == STATUS_CF_MISSING).sum()
    unnamed_n = (out_df["classification_status"] == STATUS_UNNAMED).sum()
    print(f"\n   Classification summary:")
    print(f"     fully classified (ClassyFire)     : {found}")
    print(f"     CID found, no ClassyFire data     : {no_cf}")
    print(f"     compound name not in PubChem      : {no_cid}")
    print(f"     unnamed / unidentified features   : {unnamed_n}")
    print("-- Step 2c complete ------------------------------------------------")

    return out_df


def main():
    run(config)


if __name__ == "__main__":
    main()
