"""
pipeline.py - End-to-end GCMS data processing pipeline.

Runs all processing steps in sequence:

  Step 1  data_import.py       -> output/peak_matrix_raw.csv
  Step 2  blank_correction.py  -> output/peak_matrix_blank_corrected.csv
  Step 3  normalization.py     -> output/peak_matrix_processed.csv
  Step 4  pca.py               -> output/plots/pca_scores.png + loadings
  Step 5  hca.py               -> output/plots/hca_heatmap.png + dendrogram orders

All parameters are read from config.py. Edit config.py to adapt the
pipeline to a new sample series without changing any processing code.

Usage:
    python pipeline.py

To run individual steps:
    python data_import.py
    python blank_correction.py
    python normalization.py
    python pca.py
    python hca.py
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
import normalization
import pca as pca_step
import hca as hca_step


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
    normalization.run(config)
    print()
    pca_step.run(config)
    print()
    hca_step.run(config)

    print()
    print("=" * 62)
    print("  Pipeline complete.")
    print(f"  Processed matrix : {os.path.join(config.OUTPUT_DIR, 'peak_matrix_processed.csv')}")
    print(f"  Plots            : {os.path.join(config.OUTPUT_DIR, 'plots', '')}")
    print("=" * 62)


if __name__ == "__main__":
    main()
