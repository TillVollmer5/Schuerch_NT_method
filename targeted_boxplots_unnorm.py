import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import ttest_ind

import config as cfg
import os

# --- Load blank-corrected (unnormalized) matrix ---
# Format: rows = features, columns = feature_id, S-R1..S-R6, S1..S6
matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
peak_matrix = pd.read_csv(matrix_path)

targeted_list = cfg.EXCLUSION_LIST
sample_cols_all = [c for c in peak_matrix.columns if c != "feature_id"]
column_names = ["rt", "mz", "name"] + sample_cols_all

# Build targeted data: for each target compound, find matching feature row
result_rows = []
for target in targeted_list:
    mz_target = float(target[1])
    rt_target = float(target[0])
    row_data = list(target[:3])  # rt, mz, name
    for _, feat_row in peak_matrix.iterrows():
        fid = feat_row["feature_id"].split("_")
        mz_feature = float(fid[0])
        rt_feature = float(fid[1])
        if (abs(rt_target - rt_feature) <= cfg.EXCLUSION_RT_MARGIN and
                abs(mz_target - mz_feature) <= cfg.EXCLUSION_MZ_TOLERANCE):
            row_data = list(target[:3]) + [feat_row[c] for c in sample_cols_all]
            break
    result_rows.append(row_data)

targeted_df = pd.DataFrame(result_rows, columns=column_names)
targeted_df.to_csv(os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_unnorm.csv"), index=False)


# --- Boxplots ---

sr_cols = [c for c in sample_cols_all if c.startswith("S-R")]
s_cols  = [c for c in sample_cols_all if c.startswith("S") and not c.startswith("S-R")]


def get_pvalue_signif(p_value):
    if p_value < 0.0001:
        return "****"
    elif p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"


nrows = 3
ncols = math.ceil(len(targeted_df) / nrows)
fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 2, nrows * 3))
axs = axs.flatten()
for _i in range(len(targeted_df), nrows * ncols):
    axs[_i].set_visible(False)

fig.supylabel('y-axis: unnormalized area')

for idx, row in targeted_df.iterrows():
    compound_name = row["name"]

    SR = np.array(row[sr_cols].values, dtype=float)
    S  = np.array(row[s_cols].values, dtype=float)
    print(f"{compound_name} - S count: {len(S)}, SR count: {len(SR)}")

    t_stat, p_value = ttest_ind(S, SR)
    p_signif = get_pvalue_signif(p_value)

    ax = axs[idx]

    bp = ax.boxplot([S, SR], tick_labels=['S', 'SR'], patch_artist=True)

    colors = ['coral', 'steelblue']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    for median in bp['medians']:
        median.set_color('black')

    all_data = np.concatenate([S, SR])
    if len(all_data) > 0:
        data_min = np.min(all_data)
        data_max = np.max(all_data)
        padding = (data_max - data_min) * 0.1
        ax.set_ylim(data_min - padding, data_max + padding)

    ax.set_title(f'{compound_name}\n{p_signif}')
    ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
out_path = os.path.join(cfg.OUTPUT_DIR, "plots", "targeted_boxplots_unnorm.png")
os.makedirs(os.path.join(cfg.OUTPUT_DIR, "plots"), exist_ok=True)
plt.savefig(out_path, dpi=200, bbox_inches="tight")
#plt.show()
