"""
top_features_analysis.py - Extract highest-loading PCA features with sample details.

Creates a CSV table showing the top N features by 2D Euclidean distance from
the origin in the PCA_PLOT_X / PCA_PLOT_Y loading plane — the same selection
criterion used by the loading bar chart (pca_loadings_bar.png) and loading
scatter plot (pca_loadings.png).

The number of features is controlled by PCA_BAR_TOP in config.py, so the
CSV always lists exactly the features shown in the bar chart.

Each row in the output represents one feature-sample combination where the
feature was detected (area > 0), with columns:
  - feature_id, compound_name, RT, area, sample
  - one loading column per computed PC (PC1, PC2, … up to N_COMPONENTS)

Output file: output/top_features_analysis.csv

Usage:
    python top_features_analysis.py [-n N]

    -n, --num-features    Override number of top features (default: PCA_BAR_TOP from config.py)
"""

import os
import sys
import argparse
import glob

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import pandas as pd
import numpy as np

import config


def find_compound_name(rt, peak_matrix, feature_id, data_dir, rt_tolerance=0.1):
    """
    Find compound name in raw data by matching retention time.
    Uses peak matrix to identify which samples have the feature,
    then looks up the compound name from those samples.
    
    Parameters
    ----------
    rt : float
        Retention time to match
    peak_matrix : DataFrame
        Peak matrix to identify which samples have the feature
    feature_id : str
        Feature ID to look up in peak matrix
    data_dir : str
        Directory containing raw CSV files
    rt_tolerance : float
        RT window tolerance in minutes
    
    Returns
    -------
    str or None
        Compound name if found, else None
    """
    if feature_id not in peak_matrix.index:
        return None
    
    # Get samples where this feature was detected
    feature_row = peak_matrix.loc[feature_id]
    samples_with_feature = [s for s, area in feature_row.items() if area > 0]
    
    if not samples_with_feature:
        return None
    
    # Search for compound in the samples that have this feature
    best_match = None
    best_distance = rt_tolerance
    
    for sample_name in samples_with_feature:
        # Construct expected filename
        for fp in sorted(glob.glob(os.path.join(data_dir, "*.csv"))):
            file_stem = os.path.splitext(os.path.basename(fp))[0]
            
            if file_stem != sample_name:
                continue
            
            try:
                df = None
                for skiprows in (1, 0):
                    try:
                        candidate = pd.read_csv(fp, skiprows=skiprows, engine="python", on_bad_lines="warn")
                        if "Component Name" in candidate.columns and "Retention Time" in candidate.columns:
                            df = candidate
                            break
                    except Exception:
                        pass
                if df is None:
                    continue
                
                # Find the closest RT match in this sample
                df_filtered = df[
                    (df["Retention Time"] >= rt - rt_tolerance) &
                    (df["Retention Time"] <= rt + rt_tolerance)
                ]
                
                if len(df_filtered) > 0:
                    # Find the row with closest RT
                    closest_idx = (df_filtered["Retention Time"] - rt).abs().idxmin()
                    distance = abs(df_filtered.loc[closest_idx, "Retention Time"] - rt)
                    
                    if distance < best_distance:
                        best_distance = distance
                        best_match = df_filtered.loc[closest_idx, "Component Name"]
            except Exception:
                continue
        
        if best_match:
            return best_match
    
    return None


def get_top_features(loadings_df, feature_metadata_df, n_features=10,
                     col_x="PC1", col_y="PC2"):
    """
    Get the top N features by 2D Euclidean distance from the origin in the
    col_x / col_y loading plane — identical to the selection used by the
    loading bar chart and loading scatter plot.

    Parameters
    ----------
    loadings_df : DataFrame
        Feature loadings (features x PCs).
    feature_metadata_df : DataFrame
        Columns: feature_id, mean_rt, mean_mz
    n_features : int
        Number of top features to return (config: PCA_BAR_TOP).
    col_x, col_y : str
        PC column names to use (e.g. "PC1", "PC2"); set from PCA_PLOT_X/Y.

    Returns
    -------
    DataFrame
        Top N features sorted by 2D loading distance (descending).
    """
    lx = loadings_df[col_x].values
    ly = loadings_df[col_y].values
    loadings_df['magnitude'] = np.sqrt(lx ** 2 + ly ** 2)

    # Sort by magnitude and get top N
    top = loadings_df.nlargest(n_features, 'magnitude')

    return top


def run(cfg=None, data_dir=None, output_dir=None, n_features=None):
    """
    Main analysis function.

    Parameters
    ----------
    cfg : module or None
        Config module (e.g. import config; run(config)). When provided,
        data_dir, output_dir, n_features, and PC axes are read from it.
    data_dir : str or None
        Directory containing raw sample CSV files (overrides cfg).
    output_dir : str or None
        Directory for output file (overrides cfg).
    n_features : int or None
        Number of top features to analyze. Defaults to cfg.PCA_BAR_TOP
        so the CSV always lists the same features as the bar chart.
    """
    if cfg is not None and not isinstance(cfg, str):
        # Called as run(config) from pipeline
        if data_dir is None:
            data_dir = cfg.DATA_DIR
        if output_dir is None:
            output_dir = cfg.OUTPUT_DIR
        if n_features is None:
            n_features = cfg.PCA_BAR_TOP
        pc_x = getattr(cfg, "PCA_PLOT_X", 1)
        pc_y = getattr(cfg, "PCA_PLOT_Y", 2)
    else:
        # Called standalone or as run(data_dir, output_dir)
        if isinstance(cfg, str):
            # positional: run("path/to/data")
            data_dir = cfg
        if data_dir is None:
            data_dir = config.DATA_DIR
        if output_dir is None:
            output_dir = config.OUTPUT_DIR
        if n_features is None:
            n_features = config.PCA_BAR_TOP
        pc_x = getattr(config, "PCA_PLOT_X", 1)
        pc_y = getattr(config, "PCA_PLOT_Y", 2)

    col_x = f"PC{pc_x}"
    col_y = f"PC{pc_y}"

    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print("  Top PCA Features Analysis")
    print("=" * 70)
    print(f"  Input files:")
    print(f"    {os.path.join(output_dir, 'pca_loadings.csv')}")
    print(f"    {os.path.join(output_dir, 'peak_matrix_raw.csv')}")
    print(f"    {os.path.join(output_dir, 'feature_metadata.csv')}")
    print(f"  PC axes        : {col_x} / {col_y}")
    print(f"  Number of top features: {n_features}")
    print()
    
    # Load required files
    print("  Loading data...")
    loadings_df = pd.read_csv(os.path.join(output_dir, "pca_loadings.csv"), index_col=0)
    peak_matrix = pd.read_csv(os.path.join(output_dir, "peak_matrix_raw.csv"), index_col=0)
    metadata_df = pd.read_csv(os.path.join(output_dir, "feature_metadata.csv"), index_col=0)
    
    print(f"    Loaded peak matrix with {peak_matrix.shape[0]} features x {peak_matrix.shape[1]} samples")
    print()
    
    # Get top features by loading magnitude
    print(f"  Finding top {n_features} features by loading magnitude...")
    top_features = get_top_features(loadings_df.copy(), metadata_df,
                                    n_features=n_features, col_x=col_x, col_y=col_y)
    print(f"    Found {len(top_features)} features")
    print()
    
    # Build output table
    print("  Building output table...")
    rows = []
    
    # all PC columns present in the loadings file (e.g. PC1, PC2, PC3)
    pc_cols = [c for c in loadings_df.columns if c.startswith("PC")]

    for feature_id in top_features.index:
        if feature_id not in metadata_df.index:
            continue

        mean_rt = metadata_df.loc[feature_id, "mean_rt"]
        pc_vals = {c: loadings_df.loc[feature_id, c] for c in pc_cols}

        # Find compound name
        compound_name = find_compound_name(mean_rt, peak_matrix, feature_id, data_dir)
        if compound_name is None:
            compound_name = f"Unknown (RT {mean_rt:.4f})"

        # Get areas for each sample
        if feature_id in peak_matrix.index:
            areas = peak_matrix.loc[feature_id]

            # Add one row per sample with non-zero area
            for sample_name, area in areas.items():
                if area > 0:  # Only include samples where feature was detected
                    rows.append({
                        "feature_id":    feature_id,
                        "compound_name": compound_name,
                        "RT":            mean_rt,
                        "area":          area,
                        "sample":        sample_name,
                        **pc_vals,
                    })
    
    output_df = pd.DataFrame(rows)
    
    # Calculate 2D loading distance for sorting (same metric as the plots)
    output_df['loading_magnitude'] = np.sqrt(
        output_df[col_x]**2 + output_df[col_y]**2
    )
    
    # Sort by loading magnitude (descending) then by area (descending)
    output_df = output_df.sort_values(["loading_magnitude", "area"], ascending=[False, False])
    
    # Drop the helper column
    output_df = output_df.drop(columns=['loading_magnitude'])
    
    # Write output
    output_file = os.path.join(output_dir, "top_features_analysis.csv")
    output_df.to_csv(output_file, index=False, float_format="%.6g")
    
    print(f"  ✓ Output written to: {output_file}")
    print(f"    Total rows: {len(output_df)}")
    print(f"    Unique features: {output_df['feature_id'].nunique()}")
    print()
    
    # Print summary
    print("  Top features summary:")
    print()
    for feature_id in output_df['feature_id'].unique()[:5]:
        subset = output_df[output_df['feature_id'] == feature_id]
        compound = subset.iloc[0]['compound_name']
        rt = subset.iloc[0]['RT']
        n_samples = len(subset)
        print(f"    {feature_id:12} {compound:30} RT={rt:.4f}  ({n_samples} samples)")
    
    if len(output_df['feature_id'].unique()) > 5:
        print(f"    ... and {len(output_df['feature_id'].unique()) - 5} more features")
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract highest-loading PCA features with sample details"
    )
    parser.add_argument(
        "-n", "--num-features",
        type=int,
        default=None,
        help="Number of top features to extract (default: PCA_BAR_TOP from config.py)"
    )
    
    args = parser.parse_args()
    
    run(n_features=args.num_features)
