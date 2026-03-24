"""
top_features_analysis.py - Extract highest-loading PCA features with sample details.

Creates a CSV table showing the top N features by loading magnitude, with:
  - Compound name (from raw data matching)
  - Retention time
  - Peak area in each sample
  - Sample name
  - PC1, PC2, PC3 loadings

The output includes one row per feature-sample combination where the feature
was detected (area > 0).

Output file: output/top_features_analysis.csv

Usage:
    python top_features_analysis.py [--n 10]
    
    --n, --num-features    Number of top features to extract (default: 10)
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
                df = pd.read_csv(fp, quotechar='"')
                if "Component Name" not in df.columns or "Retention Time" not in df.columns:
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


def get_top_features(loadings_df, feature_metadata_df, n_features=10):
    """
    Get the top N features by loading magnitude.
    
    Parameters
    ----------
    loadings_df : DataFrame
        Columns: PC1, PC2, PC3 (loadings for each feature)
    feature_metadata_df : DataFrame
        Columns: feature_id, mean_rt, mean_mz
    n_features : int
        Number of top features to return
    
    Returns
    -------
    DataFrame
        Top N features with feature_id, mean_rt, and loading magnitude
    """
    # Calculate magnitude of loading vector
    loadings_df['magnitude'] = np.sqrt(
        loadings_df['PC1']**2 + loadings_df['PC2']**2 + loadings_df['PC3']**2
    )
    
    # Sort by magnitude and get top N
    top = loadings_df.nlargest(n_features, 'magnitude')
    
    return top


def run(data_dir=config.DATA_DIR, output_dir=config.OUTPUT_DIR, 
        n_features=10):
    """
    Main analysis function.
    
    Parameters
    ----------
    data_dir : str
        Directory containing raw sample CSV files
    output_dir : str
        Directory for output file
    n_features : int
        Number of top features to analyze
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print("  Top PCA Features Analysis")
    print("=" * 70)
    print(f"  Input files:")
    print(f"    {os.path.join(output_dir, 'pca_loadings.csv')}")
    print(f"    {os.path.join(output_dir, 'peak_matrix_raw.csv')}")
    print(f"    {os.path.join(output_dir, 'feature_metadata.csv')}")
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
    top_features = get_top_features(loadings_df.copy(), metadata_df, n_features=n_features)
    print(f"    Found {len(top_features)} features")
    print()
    
    # Build output table
    print("  Building output table...")
    rows = []
    
    for feature_id in top_features.index:
        if feature_id not in metadata_df.index:
            continue
        
        mean_rt = metadata_df.loc[feature_id, "mean_rt"]
        pc1 = loadings_df.loc[feature_id, "PC1"]
        pc2 = loadings_df.loc[feature_id, "PC2"]
        pc3 = loadings_df.loc[feature_id, "PC3"]
        
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
                        "feature_id": feature_id,
                        "compound_name": compound_name,
                        "RT": mean_rt,
                        "area": area,
                        "sample": sample_name,
                        "PC1": pc1,
                        "PC2": pc2,
                        "PC3": pc3,
                    })
    
    output_df = pd.DataFrame(rows)
    
    # Calculate loading magnitude for sorting
    output_df['loading_magnitude'] = np.sqrt(
        output_df['PC1']**2 + output_df['PC2']**2 + output_df['PC3']**2
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
        default=10,
        help="Number of top features to extract (default: 10)"
    )
    
    args = parser.parse_args()
    
    run(n_features=args.num_features)
