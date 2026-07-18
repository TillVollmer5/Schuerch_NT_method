"""
sample_scaling.py - Step 2a of the GCMS processing pipeline.

Applies a per-sample scaling correction derived from a known physical value
(e.g. final resuspension/concentration volume in uL) supplied in config.py.
This corrects for a KNOWN technical factor, as distinct from the DATA-DRIVEN
statistical normalization methods (pqn/sum/median) applied later in
normalization.py.

Runs after blank_correction.py (so blank-vs-sample fold-change comparisons
stay on true raw instrument signal) and before normalization.py (so
pqn/sum/median can still correct for whatever technical variance the
physical value doesn't explain).

How the factor is computed
---------------------------
Raw values (e.g. "22", "5" uL) are NOT used as multipliers directly - that
would make the correction's magnitude depend on which unit was recorded in.
Instead each sample's value is expressed RELATIVE to a reference (mean or
median of all supplied values, or a fixed number):

    relative_factor = value / reference_value

SAMPLE_SCALING_DIRECTION controls how the factor is applied:
  "dilution" (default) - larger value means more dilute -> lower
                          concentration for the same underlying amount.
                          Area is MULTIPLIED by relative_factor.
  "amount"              - larger value directly means more material (e.g.
                          dry weight, cell count). Area is DIVIDED by
                          relative_factor.

Samples with no entry in SAMPLE_SCALING_VALUES get factor = 1.0 (no
correction) and a warning is printed.

Input  : output/peak_matrix_blank_corrected.csv

Output (only written when ENABLE_SAMPLE_SCALING = True):
  output/peak_matrix_blank_corrected.csv          - OVERWRITTEN with scaled
                                                     values, so every
                                                     downstream step (
                                                     normalization, boxplots,
                                                     volcano, classification,
                                                     ...) picks up the
                                                     correction transparently
  output/peak_matrix_blank_corrected_unscaled.csv - backup of the
                                                     pre-scaling matrix
  output/sample_scaling_factors.csv               - audit log: sample,
                                                     raw_value,
                                                     reference_value,
                                                     direction, factor

When ENABLE_SAMPLE_SCALING = False (default), this step is a no-op.

Usage:
    python sample_scaling.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import pandas as pd

import config


def _reference_value(values, reference):
    if isinstance(reference, (int, float)):
        return float(reference)
    series = pd.Series(values, dtype=float)
    if reference == "median":
        return float(series.median())
    if reference != "mean":
        raise ValueError(
            f"SAMPLE_SCALING_REFERENCE must be 'mean', 'median', or a number, got {reference!r}"
        )
    return float(series.mean())


def run(cfg=config):
    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 2a: sample scaling (volume/concentration correction) -----")

    enabled = getattr(cfg, "ENABLE_SAMPLE_SCALING", False)
    if not enabled:
        print("  ENABLE_SAMPLE_SCALING = False -> skipping (areas unchanged)")
        return None

    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"{matrix_path} not found - run blank_correction.py first.")

    matrix = pd.read_csv(matrix_path, index_col="feature_id")

    values_map        = getattr(cfg, "SAMPLE_SCALING_VALUES", {}) or {}
    direction          = getattr(cfg, "SAMPLE_SCALING_DIRECTION", "dilution")
    reference_setting  = getattr(cfg, "SAMPLE_SCALING_REFERENCE", "mean")

    if direction not in ("dilution", "amount"):
        raise ValueError(
            f"SAMPLE_SCALING_DIRECTION must be 'dilution' or 'amount', got {direction!r}"
        )

    known_values = [float(v) for s, v in values_map.items()
                    if s in matrix.columns and v is not None and float(v) > 0]
    if not known_values:
        print("  [warning] SAMPLE_SCALING_VALUES has no valid entries matching samples in "
              "the matrix -> skipping (areas unchanged)")
        return None

    reference_value = _reference_value(known_values, reference_setting)
    print(f"  direction     : {direction}")
    print(f"  reference     : {reference_setting}  (= {reference_value:.4g})")

    missing = [s for s in matrix.columns if s not in values_map]
    if missing:
        print(f"  [warning] no scaling value for {len(missing)} sample(s) -> factor 1.0 "
              f"(no correction): {missing}")

    audit_rows = []
    scaled = matrix.copy()
    for sample in matrix.columns:
        raw_value = values_map.get(sample)
        if raw_value is None or float(raw_value) <= 0:
            if raw_value is not None:
                print(f"  [warning] invalid scaling value for {sample!r} ({raw_value!r}) "
                      f"-> factor 1.0")
            factor = 1.0
        else:
            raw_value = float(raw_value)
            relative  = raw_value / reference_value
            factor    = relative if direction == "dilution" else (1.0 / relative)

        scaled[sample] = matrix[sample] * factor
        audit_rows.append({
            "sample":          sample,
            "raw_value":       raw_value,
            "reference_value": reference_value,
            "direction":       direction,
            "factor":          factor,
        })

    # backup pre-scaling matrix, then overwrite the canonical file so every
    # downstream step picks up the correction transparently
    backup_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected_unscaled.csv")
    matrix.to_csv(backup_path)
    print(f"  -> {backup_path}  (backup of pre-scaling matrix)")

    scaled.to_csv(matrix_path)
    print(f"  -> {matrix_path}  (overwritten with scaled values)")

    audit_df   = pd.DataFrame(audit_rows)
    audit_path = os.path.join(cfg.OUTPUT_DIR, "sample_scaling_factors.csv")
    audit_df.to_csv(audit_path, index=False)
    print(f"  -> {audit_path}")

    return scaled


if __name__ == "__main__":
    run()
