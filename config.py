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

USE_MZ       = False   # if True, also require m/z proximity to merge peaks
                       # recommended for samples with many co-eluting compounds
MZ_TOLERANCE = 0.005   # Da  - used only when USE_MZ = True



ALIGN_RT          = True   # apply median RT-shift correction across samples before detection
MZ_ALIGN_TOLERANCE = 0.1   # Da  - m/z window used during RT alignment

VALUE_COL = "Area"     # column to extract from raw CSV: "Area" or "Height"

# --- Exclusion list (normalization.py -> PCA only) ---------------------------
# Retention times (in minutes) of biologically relevant features to withhold
# from PCA so they do not dominate the principal components.
# Applied in normalization.py when building peak_matrix_processed_pca.csv.
# HCA and the volcano plot always use the full feature set (peak_matrix_processed.csv)
# so that known compounds remain visible and statistically testable there.
#
# Example entries (uncomment or add your own):
#   3.086,   # Chloroiodomethane
#   4.044,   # 3-Hexenal

EXCLUSION_LIST = [19.071,
13.587,
14.010,
26.229,
4.043,
5.643,
5.690,
10.614,
23.321,
22.810,
8.090,
9.529,
11.290,
10.591,
11.918,
10.072

    # 3.086,
]

EXCLUSION_RT_MARGIN = 0.05   # +- minutes around each listed RT

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
MIN_PREVALENCE_HCA     = 0.0   # set > 0 to drop sparse features from the heatmap
MIN_PREVALENCE_VOLCANO = 0.0   # leave at 0.0 to keep group-specific features

# --- Blank correction (blank_correction.py) -----------------------------------
FOLD_CHANGE_THRESHOLD = 10.0
# A feature is removed if:
#   mean(sample areas) / max(blank area)  <  FOLD_CHANGE_THRESHOLD
# Features absent from blanks are always retained.
# Common values: 3 (lenient) - 10 (strict)

# --- Normalization (normalization.py) -----------------------------------------
NORMALIZATION = "sum"
# "sum"    - divide each sample by its total signal (scaled to median column sum)
# "median" - scale each sample so its median equals the global median
# "none"   - skip this step

LOG_BASE = 2
# Logarithm base for transformation: 2 (log2) or math.e (natural log)
# log2 is conventional in metabolomics

SCALING = "pareto"
# "pareto" - mean-centre, divide by sqrtstd  (recommended for untargeted GCMS)
# "auto"   - mean-centre, divide by std   (unit-variance; stronger equalisation)
# "none"   - skip scaling (apply only after log transform)

# --- Volcano plot (volcano.py) -----------------------------------------------
VOLCANO_COMPARISONS  = "all"
# "all" - run every pairwise group comparison automatically
# or specify a list of (groupA, groupB) tuples, e.g.:
#   VOLCANO_COMPARISONS = [("S", "S-R")]
# log2FC is expressed as groupA / groupB (positive = higher in A)

VOLCANO_FC_THRESHOLD  = 1.0    # log2 fold-change cutoff (1.0 = 2-fold change)
VOLCANO_P_THRESHOLD   = 0.05   # Benjamini-Hochberg adjusted p-value threshold
VOLCANO_TOP_LABELS    = 10     # number of top significant features to label in the plot
                                # ranked by adjusted p-value; set to 0 to suppress labels

# --- HCA (hca.py) -------------------------------------------------------------
HCA_LINKAGE           = "ward"       # linkage method: "ward", "average", "complete", "single"
                                     # note: ward requires metric="euclidean"
HCA_METRIC            = "euclidean"  # distance metric; switch to "correlation" with
                                     # linkage "average" or "complete"
HCA_CMAP              = "vlag"       # diverging colormap suited to mean-centred scaled data
HCA_MAX_FEATURE_LABELS = 50          # label the feature axis when n_features <= this value;
                                     # set to 0 to always hide feature labels

# --- PCA (pca.py) -------------------------------------------------------------
N_COMPONENTS    = 3   # number of principal components to compute and save
                      # increase to retain more dimensions (e.g. 5 for scree plot)

PCA_PLOT_X      = 1   # PC number to plot on the X axis (1-indexed)
PCA_PLOT_Y      = 2   # PC number to plot on the Y axis (1-indexed)

PCA_ELLIPSE      = True  # draw 95 % confidence ellipses per group (requires scipy)
PCA_TOP_LOADINGS = 10   # number of top-loading features to label in the loadings scatter plot
                        # set to 0 to skip labels
PCA_BAR_TOP      = 10   # number of features shown in the loading bar chart (pca_loadings_bar.png)
                        # selected by Euclidean distance in the PC_x/PC_y loading plane;
                        # increase to inspect more candidates (e.g. 20)


