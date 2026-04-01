"""
pipeline.py - End-to-end GCMS data processing pipeline.

Runs all processing steps in sequence:

  Step 1  data_import.py            -> output/peak_matrix_raw.csv
  Step 2  blank_correction.py       -> output/peak_matrix_blank_corrected.csv
  Step 2b prevalence_histogram.py   -> output/plots/prevalence_histogram.png
  Step 2c compound_classification.py -> output/compound_classes.csv (optional,
                                        controlled by RUN_COMPOUND_CLASSIFICATION)
  Step 2d compound_class_plots.py   -> output/plots/class_pie_*.png (optional,
                                        controlled by RUN_CLASS_PLOTS)
  Step 3  normalization.py          -> output/peak_matrix_processed.csv
  Step 4  pca.py                    -> output/plots/pca_scores.png + loadings
  Step 5  hca.py                    -> output/plots/hca_heatmap.png + dendrogram orders
  Step 5b hca_dendrogram.py         -> output/plots/hca_dendrogram.html (interactive)
  Step 6  volcano.py                -> output/plots/volcano_*.png + results tables
  Step 7  top_features_analysis.py  -> output/top_features_analysis.csv
  Step 8  blank_contaminants_report.py -> output/blank_contaminants_report.csv

All parameters are read from config.py. Edit config.py to adapt the
pipeline to a new sample series without changing any processing code.

Usage:
    python pipeline.py

To run individual steps:
    python data_import.py
    python blank_correction.py
    python prevalence_histogram.py
    python compound_classification.py   # optional; respects RUN_COMPOUND_CLASSIFICATION
    python compound_class_plots.py      # optional; respects RUN_CLASS_PLOTS
    python normalization.py
    python pca.py
    python hca.py
    python volcano.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import config
import data_import
import blank_correction
import prevalence_histogram as prevalence_histogram_step
import compound_classification as compound_classification_step
import compound_class_plots as compound_class_plots_step
import normalization
import pca as pca_step
import hca as hca_step
import hca_dendrogram as hca_dendrogram_step
import volcano as volcano_step
import top_features_analysis as top_features_analysis_step
import blank_contaminants_report as blank_contaminants_report_step


def main():
    print("=" * 62)
    print("  GCMS data processing pipeline")
    print(f"  data dir     : {config.DATA_DIR}")
    print(f"  output dir   : {config.OUTPUT_DIR}")
    print(f"  RT margin    : {config.RT_MARGIN} min")
    print(f"  use m/z      : {config.USE_MZ}")
    print(f"  fold change  : {config.FOLD_CHANGE_THRESHOLD}x")
    print(f"  normalization: {config.NORMALIZATION}")
    print(f"  log base     : {config.LOG_BASE}")
    print(f"  scaling      : {config.SCALING}")
    print(f"  PCA components: {config.N_COMPONENTS}")
    print("=" * 62)

    data_import.run(config)
    print()
    blank_correction.run(config)
    print()
    prevalence_histogram_step.run(config)
    print()
    if getattr(config, "RUN_COMPOUND_CLASSIFICATION", False):
        compound_classification_step.run(config)
        print()
    if getattr(config, "RUN_CLASS_PLOTS", False):
        compound_class_plots_step.run(config)
        print()
    normalization.run(config)
    print()
    pca_step.run(config)
    print()
    hca_step.run(config)
    print()
    if getattr(config, "RUN_HCA_DENDROGRAM", True):
        hca_dendrogram_step.run(config)
        print()
    volcano_step.run(config)
    print()
    top_features_analysis_step.run(config)
    print()
    blank_contaminants_report_step.run(config)

    print()
    print("=" * 62)
    print("  Pipeline complete.")
    print(f"  Processed matrix : {os.path.join(config.OUTPUT_DIR, 'peak_matrix_processed.csv')}")
    print(f"  Plots            : {os.path.join(config.OUTPUT_DIR, 'plots', '')}")
    print("=" * 62)


if __name__ == "__main__":
    main()
