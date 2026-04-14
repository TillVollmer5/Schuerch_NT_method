"""
config.py - Central configuration for the GCMS processing pipeline.

Edit this file to adapt the pipeline to a new sample series.
All scripts import their parameters from here.
"""

import math

# --- File paths ---------------------------------------------------------------
DATA_DIR   = "DATA"     # directory containing raw TraceFinder CSV exports
OUTPUT_DIR = "output"   # all output files are written here

# --- Sample groups ------------------------------------------------------------
# Files whose names start with BLANK_PREFIX are treated as blanks (excluded
# from the main matrix but used for blank correction).
#
# All other files are classified into groups by matching their filename prefix
# against SAMPLE_GROUPS in order - put more specific prefixes first so that
# e.g. "S-R" is matched before "S".
#
# To adapt to a new sample series, change the group names and/or prefixes here.
# Example for a three-group experiment:
#   SAMPLE_GROUPS = [("control", "C"), ("treatment_A", "TA"), ("treatment_B", "TB")]

BLANK_PREFIX  = "Blank"

SAMPLE_GROUPS = [
    ("S-R", "S-R"),   # group 1 - reference/treatment  (S-R1, S-R2, S-R3)
    ("S",   "S"),     # group 2 - samples               (S1-S6)
]

# --- Feature detection (data_import.py) --------------------------------------
RT_MARGIN    = 0.05    # minutes  - peaks within this RT window -> same feature
                       # increase for lower RT reproducibility, decrease for dense spectra

USE_MZ       = True   # if True, also require m/z proximity to merge peaks
                       # recommended for samples with many co-eluting compounds
MZ_TOLERANCE = 0.0005   # Da  - used only when USE_MZ = True

ALIGN_RT          = True   # apply median RT-shift correction across samples before detection
MZ_ALIGN_TOLERANCE = 0.1   # Da  - m/z window used during RT alignment

VALUE_COL = "Area"     # column to extract from raw CSV: "Area" or "Height"

# --- Exclusion list (normalization.py -> PCA only) = Targeted list ---------------------------
# Retention times (in minutes) and m/z values of biologically relevant features to withhold
# from PCA so they do not dominate the principal components.
# Applied in normalization.py when building peak_matrix_processed_pca.csv.
# HCA and the volcano plot always use the full feature set (peak_matrix_processed.csv)
# so that known compounds remain visible and statistically testable there.
#
# Example entries (uncomment or add your own):
#   [5.997, 41.0384, "Z-3-Hexenal"],#3-Hexenal
#   [7.735, 83.0492, "E-2-Hexenal"],#2-Hexenal, (E)-
#   [rt,mz,"name"]  - specify RT, m/z, and a name for reference; m/z can be None to ignore m/z in matching

EXCLUSION_LIST = [
#[5.997, 41.0384, "Z-3-Hexenal"],#3-Hexenal
#[7.735, 83.0492, "E-2-Hexenal"],#2-Hexenal, (E)- E/Z based on literature
#[7.801, 67.0542, "Z-3-Hexenol"],#3-Hexen-1-ol, (Z)-
#[8.958, 104.0621, "Styrene"],#Styrene
#[10.43, 91.0542, "alpha-Thujene"],#.alpha.-Thujene or pinene
#[11.934,105.0699, "sesquiterpene"],#2,3-Diazabicyclo   pinene or thujene
#[12.88, 67.0542, "Z-3-Hexenol acetate"],#3-Hexen-1-ol, acetate, (Z)- e/z based on hexenal
#[13.205,67.0542, "E-2-Hexenol acetate"],#2-Hexen-1-ol, acetate, (E)-
#[13.522,119.0856, "Cymene"],#p-Cymene or any other cymene isomer
#[13.665,41.0384, "Limonene"],#limonene reference
#[14.254,91.0542, "E-beta-Ocimene"],#trans-.beta.-Ocimene
#[15.93, 93.0699, "Linalool"],#Linalool
#[16.375,41.0384, "DMNT"],#4,8-DIMETHYLNONA-1,3,7-TRIENE
#[21.529,117.0573, "Indole"],#Indole
#[25.049,91.0542, "(-)-(E)-Caryophyllene"],#(-)-(E)-Caryophyllene reference
#[25.38, 119.0856, "E-alpha-Bergamotene"],#trans-.alpha.-Bergamotene
#[25.828,91.0542, "Isogermacrene D"],#Isogermacrene D
#[28.763, 81.0699, "TMTT"]#(3E,7E)-4,8,12-Trimethyltrideca-1,3,7,11-tetraene
]

EXCLUSION_RT_MARGIN = 0.05   # +- minutes around each listed RT
EXCLUSION_MZ_TOLERANCE = MZ_TOLERANCE   # +- Da around each listed m/z based on maximal variance observed

# --- Missingness / prevalence filter ------------------------------------------
# Features coded as 0 after blank correction are "not detected" in that sample.
# These parameters control the minimum fraction of samples in which a feature
# must be detected (area > 0) to be retained in each analysis.
# Range: 0.0 (keep all) to 1.0 (require detection in every sample).
#
# PCA  - strict filter recommended: noise features with many zeros destabilise
#        principal components without adding biological information.
# HCA  - moderate filter; sparse features create uninformative columns but
#        group-shared features should be retained. Set to 0.0 to disable.
# Volcano - keep at 0.0 by default: a compound present in ALL samples of one
#        group but NONE of the other is the most biologically interesting result.
#        Filtering by overall prevalence would remove exactly those features.

MIN_PREVALENCE_PCA     = 0.35   # e.g. 0.5 = detected in >= 50% of all samples
MIN_PREVALENCE_HCA     = 0.45   # set > 0 to drop sparse features from the heatmap
MIN_PREVALENCE_VOLCANO = 0.0   # leave at 0.0 to keep group-specific features

# --- Blank correction (blank_correction.py) -----------------------------------
FOLD_CHANGE_THRESHOLD = 3.0
# A feature is removed if:
#   mean(sample areas) / max(blank area)  <  FOLD_CHANGE_THRESHOLD
# Features absent from blanks are always retained.
# Common values: 3 (lenient) - 10 (strict)

BLANK_USE_MZ       = True   # if True, also require m/z proximity for a blank peak
                              # to count as a match against a sample feature.
                              # Without this, a blank peak at the same RT but a
                              # different m/z (a different compound) can cause a
                              # sample feature to be incorrectly removed.
BLANK_MZ_TOLERANCE = MZ_TOLERANCE   # Da  - maximum |feature_mz - blank_mz| to accept match
                              # used only when BLANK_USE_MZ = True

BLANK_SAMPLE_MODE = "per_sample"
# How the sample signal is represented for the fold-change comparison:
# "mean"       - arithmetic mean across all samples (current/default behaviour)
# "per_sample" - each sample is compared individually; cells that fail are
#                zeroed in the matrix rather than removing the whole feature row.
#                A feature is fully dropped only if every sample cell is zeroed.
# "per_group"  - mean area per group (from SAMPLE_GROUPS); failing groups have
#                all their sample cells zeroed. Same all-zero -> drop rule applies.

BLANK_REFERENCE_MODE = "max"
# Which blank signal is used as the fold-change denominator:
# "max"  - highest area across all blank files within the RT window (current/default)
# "mean" - mean area across all blank files (m/z gate applied per blank before averaging)
# "each" - each blank file is compared independently; a sample/group/mean fails
#           if its fold change falls below FOLD_CHANGE_THRESHOLD for ANY blank file

BLANK_EXCLUDE_KEYWORDS = ["silan", "Silan", "Si", "siloxane", "Siloxane"]
# Features whose compound name or molecular formula contains any of these
# substrings (case-insensitive) are removed after blank correction.
# Useful for stripping known instrument/column contaminants by name or element.
# Examples:
#   BLANK_EXCLUDE_KEYWORDS = ["silan", "Si"]   # removes siloxanes / Si-containing compounds
#   BLANK_EXCLUDE_KEYWORDS = ["column", "phthalate"]

# --- Normalization, log transform, and scaling --------------------------------
#
# Global defaults — used by any analysis that does not have an override below.
# Each analysis (PCA, HCA, Volcano, Boxplot) can override these independently
# via the NORMALIZATION_<X>, LOG_BASE_<X>, and SCALING_<X> keys.
# Set an override to None to fall back to the corresponding global default.
#
# Normalization options (applied sample-wise on RAW areas, before log):
#   "none"   - skip normalization
#   "sum"    - divide each sample by its total signal (rescaled to median sum)
#   "median" - scale each sample so its median equals the global median
#   "pqn"    - Probabilistic Quotient Normalization (Dieterle et al. 2006)
#              Recommended for VOC data where total emission differs biologically.
#              Computes per-sample quotients vs a reference spectrum (feature-wise
#              median), takes the median quotient as the normalization factor.
#
# Log transform options (applied after normalization):
#   2        - log2(x+1)   conventional in metabolomics (default)
#   10       - log10(x+1)  useful for very wide dynamic ranges
#   math.e   - ln(x+1)     natural log
#   "sqrt"   - sqrt(x)     cube-root-like compression, handles zeros
#   "cbrt"   - cbrt(x)     cube-root transform, symmetric for negative values
#   "none"   - skip log transform
#
# Scaling options (applied feature-wise after log, per feature across all samples):
#   "none"   - skip scaling
#   "pareto" - (x - mean) / sqrt(std)  recommended for untargeted GC-MS
#   "auto"   - (x - mean) / std        unit-variance; stronger equalisation
#   "vast"   - (x - mean) * mean / std²  VAST: weights stable features (Van den Berg 2006)
#   "range"  - (x - mean) / (max - min)  range scaling
#   "level"  - (x - mean) / mean         level scaling (relative to mean)

NORMALIZATION = "pqn"
LOG_BASE      = 2
SCALING       = "pareto"

# --- Per-analysis overrides ---------------------------------------------------
# Set any value to None to inherit the corresponding global setting above.
# Recommended starting point for VOC / GC-MS metabolomics:
#
#   PCA  / HCA  : PQN + log2 + pareto  (pattern discovery, robust to abundance bias)
#   Volcano     : PQN + log2 + none    (IMPORTANT: keep scaling "none" here — scaling
#                                       distorts fold-change magnitudes)
#   Boxplot     : none + none + none   (raw areas; directional biology is clear)

NORMALIZATION_PCA     = None   # inherits NORMALIZATION
NORMALIZATION_HCA     = None   # inherits NORMALIZATION
NORMALIZATION_VOLCANO = None   # inherits NORMALIZATION
NORMALIZATION_BOXPLOT = "none" # raw areas recommended for targeted inspection

LOG_BASE_PCA     = None   # inherits LOG_BASE
LOG_BASE_HCA     = None   # inherits LOG_BASE
LOG_BASE_VOLCANO = None   # inherits LOG_BASE
LOG_BASE_BOXPLOT = "none" # no log for boxplots (raw areas on y-axis)

SCALING_PCA     = None    # inherits SCALING
SCALING_HCA     = None    # inherits SCALING
SCALING_VOLCANO = "none"  # must stay "none" for valid fold-change computation
SCALING_BOXPLOT = "none"  # no scaling for boxplots

# --- Targeted compound list (targeted_boxplots.py) ---------------------------
# Compounds to visualise as individual boxplots, one panel per compound.
# Independent from EXCLUSION_LIST — you can include compounds from there or add others.
#
# Each entry is a list with the following fields:
#   [rt, mz, "name"]
#       Match by retention time (+-EXCLUSION_RT_MARGIN) and m/z (+-EXCLUSION_MZ_TOLERANCE).
#
#   [rt, mz, "name", "feature_id"]
#       Match by exact feature_id string (e.g. "41.0384_5.997").
#       rt and mz are still used for display; set them to None if unknown.
#
# If TARGETED_LIST is empty ([]), the script falls back to EXCLUSION_LIST.

TARGETED_LIST = [
[5.997,  41.0384, "Z-3-Hexenal"],
[7.735,  83.0492, "E-2-Hexenal"],
[7.801,  67.0542, "Z-3-Hexenol"],
[8.958, 104.0621, "Styrene"],
[12.88,  67.0542, "Z-3-Hexenol acetate"],
[13.205, 67.0542, "E-2-Hexenol acetate"],
[13.522,119.0856, "Cymene"],
[13.665, 41.0384, "Limonene"],
[15.93,  93.0699, "Linalool"],
[16.375, 41.0384, "DMNT"],
[21.529, 117.0573, "Indole"],
[25.049, 91.0542, "(-)-(E)-Caryophyllene"],
[28.763, 81.0699, "TMTT"],
]

TARGETED_BOXPLOT_ROWS = 3
# Number of rows in the boxplot grid.  Columns are computed automatically.
# Increase for more compounds to keep panels readable.

RUN_TARGETED_BOXPLOTS = True
# Set to False to skip targeted boxplot generation in pipeline.py.

STAT_TEST_BOXPLOT = "mannwhitney"
# Statistical test for pairwise group comparison in targeted boxplots.
# Applied between the first two groups in SAMPLE_GROUPS.
#
# "mannwhitney"  - Mann-Whitney U test (non-parametric; recommended for n < 10)
#                  Tests whether one group tends to have higher values than the other.
#                  Robust to non-normality and outliers.
# "ttest"        - Welch's t-test (parametric; assumes approximate normality)
#                  Use when n ≥ 8 and data distribution is roughly normal.
# "ttest_equal"  - Student's t-test (equal variances assumed)
#                  Use only when you have reason to believe variances are equal.
# "wilcoxon"     - Wilcoxon signed-rank test (paired non-parametric)
#                  Use only when samples are paired (matched treatment/control).
# "spearman"     - Spearman rank correlation between group index and values
#                  Useful when groups represent an ordered gradient (dose series).
# "pearson"      - Pearson correlation between group index and values
#                  Use only when relationship is expected to be linear and normal.
# "kruskal"      - Kruskal-Wallis test (non-parametric ANOVA; >2 groups)
#                  Recommended when comparing three or more groups.
# "anova"        - One-way ANOVA (parametric; >2 groups)
#                  Use when normality and equal variances can be assumed.

# --- Volcano plot (volcano.py) -----------------------------------------------
VOLCANO_COMPARISONS  = "all"
# "all" - run every pairwise group comparison automatically
# or specify a list of (groupA, groupB) tuples, e.g.:
#   VOLCANO_COMPARISONS = [("S", "S-R")]
# log2FC is expressed as groupA / groupB (positive = higher in A)

VOLCANO_FC_THRESHOLD  = 1.0    # log2 fold-change cutoff (1.0 = 2-fold change)
VOLCANO_P_THRESHOLD   = 0.05   # Benjamini-Hochberg adjusted p-value threshold
VOLCANO_TOP_LABELS    = 30     # number of top significant features to label in the plot
                                # ranked by adjusted p-value; set to 0 to suppress labels

STAT_TEST_VOLCANO = "mannwhitney"
# Statistical test used to compute per-feature p-values in the volcano plot.
# The log2 fold-change (mean difference on the log scale) is always used for
# the x-axis regardless of test choice.  All p-values are BH-corrected.
# IMPORTANT: do NOT use correlation tests here (spearman/pearson) — they are
# not appropriate for two-group volcano analysis.
#
# "mannwhitney"  - Mann-Whitney U  (non-parametric; recommended default)
# "ttest"        - Welch's t-test  (parametric; larger n, normal distribution)
# "ttest_equal"  - Student's t-test (equal variances)
# "kruskal"      - Kruskal-Wallis  (extends to >2 groups if needed)

# --- HCA (hca.py) -------------------------------------------------------------
HCA_LINKAGE           = "ward"       # linkage method: "ward", "average", "complete", "single"
                                     # note: ward requires metric="euclidean"
HCA_METRIC            = "euclidean"  # distance metric; switch to "correlation" with
                                     # linkage "average" or "complete"
HCA_CMAP              = "vlag"       # diverging colormap suited to mean-centred scaled data
HCA_MAX_FEATURE_LABELS = 50          # label the feature axis when n_features <= this value;
                                     # set to 0 to always hide feature labels

RUN_HCA_DENDROGRAM = True
# Set to False to skip the interactive HTML dendrogram (Step 5b) in pipeline.py.
# Produces hca_dendrogram.html — no server needed, opens in any browser.

HCA_CLASS_ANNOTATION_COLUMNS = ["superclass", "npclassifier_pathway"]
# Columns from feature_metadata_enriched.csv to display as colored annotation
# strips alongside the feature (column) dendrogram in the HCA heatmap.
# Each entry produces one strip; strips reorder automatically with the dendrogram
# so you can see which compound classes cluster together.
# Set to [] to show no annotation strips.
#
# Recommended options for your data:
#   ["superclass"]
#   ["superclass", "npclassifier_pathway"]
#   ["superclass", "npclassifier_pathway", "subclass"]
#
# Available columns (from feature_metadata_enriched.csv):
#   "superclass"              - broad chemical class (Lipids, Benzenoids, ...)
#   "npclassifier_pathway"    - biosynthetic pathway (Terpenoids, Fatty acids, ...)
#   "subclass"                - more specific grouping (Sesquiterpenoids, ...)
#   "npclassifier_class"      - NPClassifier class level
#   "kingdom", "class", "direct_parent"  - other ClassyFire levels

# --- Feature labelling (pca.py, hca.py, volcano.py) --------------------------
COMPOUND_NAME_COL = "Component Name"   # column in TraceFinder CSVs holding compound names
                              # set to "" to disable name extraction

FEATURE_LABEL = "id"         # how to label features in plots and axes
                              # "id"   - use the auto-generated feature_id  (default)
                              # "name" - use the compound name of the most abundant
                              #          signal in the cluster; falls back to feature_id
                              #          when no name is available

# --- PCA (pca.py) -------------------------------------------------------------
N_COMPONENTS    = 2   # number of principal components to compute and save
                      # increase to retain more dimensions (e.g. 5 for scree plot)

PCA_PLOT_X      = 1   # PC number to plot on the X axis (1-indexed)
PCA_PLOT_Y      = 2   # PC number to plot on the Y axis (1-indexed)

PCA_ELLIPSE      = True  # draw 95 % confidence ellipses per group (requires scipy)
PCA_TOP_LOADINGS = 10   # number of top-loading features to label in the loadings scatter plot
                        # set to 0 to skip labels
PCA_BAR_TOP      = 10   # number of features shown in the loading bar chart (pca_loadings_bar.png)
                        # and exported to top_features_analysis.csv.
                        # Selected by Euclidean distance in the PCA_PLOT_X/PCA_PLOT_Y loading plane.

# --- Compound class plots (compound_class_plots.py) --------------------------
RUN_CLASS_PLOTS = True
# Set to False to skip pie chart generation entirely in pipeline.py.

CLASS_PIE_COLUMNS = ["superclass", "npclassifier_pathway", "npclassifier_superclass", "subclass"]
# List of columns from feature_metadata_enriched.csv to visualise as pie charts.
# Any column in that file is valid, e.g.:
#   ["kingdom", "superclass", "class", "subclass", "direct_parent",
#    "npclassifier_pathway", "npclassifier_superclass", "npclassifier_class"]

CLASS_PIE_GROUPS = "separate"
# Which samples to include in each pie:
#   "separate" - one pie per group defined in SAMPLE_GROUPS
#   "combined" - one pie across all samples
#   ["S", "S-R"] - explicit list; only the listed groups are plotted

CLASS_PIE_MIN_FRACTION = 0.02
# Slices that represent less than this fraction of total features are merged
# into an "Other" slice.  Set to 0.0 to show all slices.

CLASS_PIE_DETECTED_ONLY = True
# True  - count only features detected (area > 0) in at least one sample of
#         the group (group-specific feature presence)
# False - count all features in the enriched metadata regardless of detection

# --- Class highlighting in PCA / volcano plots --------------------------------
# Each entry colors all features matching the given column/value combination.
# Highlighted features appear as colored dots in the PCA loadings scatter,
# colored y-tick labels in the loadings bar chart, and colored ring outlines
# on significant dots in the volcano plot.
# Features matching multiple entries use the last matching color.
#
# Available superclass values in your data:
#   "Lipids and lipid-like molecules"  "Phenylpropanoids and polyketides"
#   "Benzenoids"                       "Organohalogen compounds"
#   "Hydrocarbon derivatives"          "Hydrocarbons"
#   "Organic oxygen compounds"         "Organic acids and derivatives"
#   "Organoheterocyclic compounds"     "Organosulfur compounds"
#   "Alkaloids and derivatives"        "Organophosphorus compounds"
#
# Available npclassifier_pathway values in your data:
#   "Terpenoids"   "Fatty acids"   "Polyketides"
#   "Shikimates and Phenylpropanoids"   "Alkaloids"
#
# Available subclass values (selection):
#   "Sesquiterpenoids"  "Monoterpenoids"  "Diterpenoids"  "Triterpenoids"
#   "Fatty acids and conjugates"  "Fatty alcohols"  "Fatty acid esters"
#   "Organohalogen compounds"

CLASS_HIGHLIGHT = [
    # Uncomment and edit to activate highlighting.  Use any column from
    # feature_metadata_enriched.csv and any value listed above.
    {"column": "npclassifier_pathway", "value": "Terpenoids",                          "color": "#edaf29"},
    {"column": "subclass",             "value": "Sesquiterpenoids",                    "color": "#2ecc71"},
    {"column": "subclass",             "value": "Monoterpenoids",                      "color": "#3498db"},
    {"column": "subclass",             "value": "Diterpenoids",                        "color": "#9b59b6"},
    {"column": "subclass",             "value": "Triterpenoids",                        "color": "#e67e22"},
    {"column": "subclass",             "value": "Fatty acids and conjugates",          "color": "#1abc9c"},
    {"column": "subclass",             "value": "Fatty alcohols",                        "color": "#16a085"},
    {"column": "subclass",             "value": "Fatty acid esters",                        "color": "#27ae60"},
    {"column": "subclass",             "value": "Oxanes",                        "color": "#c0392b"},
    # {"column": "superclass",           "value": "Organohalogen compounds",             "color": "#e74c3c"},
    # {"column": "superclass",           "value": "Phenylpropanoids and polyketides",    "color": "#9b59b6"},
]

CLASS_LABEL_COLUMN = "subclass"
# When non-empty, the value from this column is appended to feature labels in
# the PCA loadings scatter, loadings bar chart, and volcano plot as "[class]",
# e.g. "F042 [Sesquiterpenoids]".
# Only applied when the column has a non-NaN value for that feature.
#
# Recommended options for your data:
#   CLASS_LABEL_COLUMN = "subclass"            # most specific useful level
#   CLASS_LABEL_COLUMN = "superclass"          # broader grouping
#   CLASS_LABEL_COLUMN = "npclassifier_class"  # NP biosynthetic class

# --- Compound class color palette -------------------------------------------
# Global color registry: maps class value strings → hex colors.
# Used consistently in:
#   - Compound class pie charts    (Step 2d  compound_class_plots.py)
#   - HCA annotation strips        (Step 5   hca.py)
#   - Interactive dendrogram pies  (Step 5b  hca_dendrogram.py)
#
# Values not listed here are auto-assigned from a tab20 palette in sorted order,
# so the same unrecognised value always gets the same auto-assigned color.
# "Unknown", "Unclassified", and "Other" always receive their gray entries below.
#
# To add or override a color, add:
#   "YourClassValue": "#rrggbb",

CLASS_COLORS = {
    # ---- ClassyFire superclass -----------------------------------------------
    "Lipids and lipid-like molecules":       "#ff7f0e",
    "Phenylpropanoids and polyketides":      "#9467bd",
    "Benzenoids":                            "#1f77b4",
    "Organohalogen compounds":               "#d62728",
    "Hydrocarbon derivatives":               "#8c564b",
    "Hydrocarbons":                          "#7f7f7f",
    "Organic oxygen compounds":              "#17becf",
    "Organic acids and derivatives":         "#bcbd22",
    "Organoheterocyclic compounds":          "#e377c2",
    "Organosulfur compounds":                "#aec7e8",
    "Alkaloids and derivatives":             "#98df8a",
    "Organophosphorus compounds":            "#ffbb78",
    # ---- NPClassifier pathway ------------------------------------------------
    "Terpenoids":                            "#2ca02c",
    "Fatty acids":                           "#ff9896",
    "Polyketides":                           "#c5b0d5",
    "Shikimates and Phenylpropanoids":       "#f7b6d2",
    "Alkaloids":                             "#dbdb8d",
    # ---- ClassyFire subclass (selection) ------------------------------------
    "Sesquiterpenoids":                      "#2ecc71",
    "Monoterpenoids":                        "#3498db",
    "Diterpenoids":                          "#1abc9c",
    "Triterpenoids":                         "#27ae60",
    "Fatty acids and conjugates":            "#e74c3c",
    "Fatty alcohols":                        "#c0392b",
    "Fatty acid esters":                     "#f39c12",
    # ---- Special / fallback (always gray) -----------------------------------
    "Unknown":                               "#cccccc",
    "Unclassified":                          "#aaaaaa",
    "Other":                                 "#bbbbbb",
}

# --- Compound classification (compound_classification.py) --------------------
RUN_COMPOUND_CLASSIFICATION = True
# Set to False to skip the PubChem classification step entirely in pipeline.py.

PUBCHEM_CACHE_ONLY = True
# True  - build the output from the local cache only; no network requests are
#         made.  Compounds not yet in the cache are marked "unnamed" in the
#         output instead of being queried.  Use this when you are offline, want
#         a fast re-run, or simply want to regenerate the output CSVs without
#         consuming API quota.
# False - fetch missing entries from PubChem as normal (first run or after
#         deleting the cache file).  Already-cached entries are still served
#         from the cache and never re-fetched.

PUBCHEM_USER_AGENT = "Schuerch_NT_pipeline/1.0 (nontargeted GCMS metabolomics; contact: till.vollmer@unibe.ch)"
# Replace YOUR_EMAIL_HERE with your real email address.
# PubChem's usage policy requests a descriptive User-Agent so they can contact
# you if your script causes unexpected server load.
# Policy: https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access

PUBCHEM_RATE_LIMIT_DELAY = 1
# Seconds to sleep before each PubChem API request.
# PubChem's stated limit is 5 requests/sec; 0.35 s gives ~2.9/sec (conservative).
# Increase (e.g. to 0.5) if you receive frequent 503 responses.

CLASSYFIRE_RATE_LIMIT_DELAY = 3
# Seconds to sleep before each ClassyFire API request (classyfire.wishartlab.com).
# No stated rate limit; 2.0 s is conservative — 1.0 s has been observed to
# trigger HTTP 429 (Too Many Requests) responses from this academic server.

NPCLASSIFIER_RATE_LIMIT_DELAY = 3
# Seconds to sleep before each NPClassifier API request (npclassifier.gnps2.org).
# No stated rate limit; 1.0 s is conservative and well-mannered.

PUBCHEM_CACHE_FILE = "output/pubchem_cache.json"
# Path to the local JSON cache for API responses.
# The cache stores: compound name -> CID, CID -> properties,
#                   CID -> SMILES + InChIKey, CID -> ClassyFire taxonomy,
#                   CID -> NPClassifier annotations.
# Delete this file (or a specific entry inside it) to force a re-fetch.
