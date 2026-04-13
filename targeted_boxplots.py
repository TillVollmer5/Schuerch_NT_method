import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import ttest_ind

import config as cfg
import os
matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed.csv")
peak_matrix_processed = pd.read_csv(matrix_path)
targeted_list=cfg.EXCLUSION_LIST
column_names=["rt","mz","name"]

for row in peak_matrix_processed["sample"]:
    column_names.append(row)




# Iterate over column names of full processed data matrix
iter=0 #first column is samples names, so skip it
for column in peak_matrix_processed:
    if iter>0:
        #see if feature ids match
        feature_id=column.split("_")
        mz_feature=float(feature_id[0])
        rt_feature=float(feature_id[1])
        for row in targeted_list:
            mz_target=float(row[1])
            rt_target=float(row[0])
            if abs(rt_target - rt_feature) <= cfg.EXCLUSION_RT_MARGIN and abs(mz_target - mz_feature) <= cfg.EXCLUSION_MZ_TOLERANCE:
                column_list=peak_matrix_processed[column].tolist()
                for i in range(len(column_list)):
                    row.append(column_list[i])
    iter+=1

targeted_df=pd.DataFrame(targeted_list)
targeted_df.rename(columns={i: column_names[i] for i in range(len(column_names))}, inplace=True)
pd.DataFrame(targeted_df).to_csv(os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot.csv"), index=False)



# --- Boxplots ---

# Identify sample columns by group from column names
sample_cols = [c for c in targeted_df.columns if c not in ("rt", "mz", "name")]
sr_cols = [c for c in sample_cols if c.startswith("S-R")]
s_cols  = [c for c in sample_cols if c.startswith("S") and not c.startswith("S-R")]

def get_pvalue_signif(p_value):
    # Convert p-value to significance symbols (like R's p.signif)
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
# hide unused subplots if grid is larger than number of compounds
for _i in range(len(targeted_df), nrows * ncols):
    axs[_i].set_visible(False)

fig.supylabel('y-axis: normalized area')

# Create boxplots for each compound
for idx, row in targeted_df.iterrows():
    compound_name = row["name"]

    SR = np.array(row[sr_cols].values, dtype=float)
    S  = np.array(row[s_cols].values, dtype=float)
    print(f"{compound_name} - S count: {len(S)}, SR count: {len(SR)}")

    # Perform two-sample t-test and get significance
    t_stat, p_value = ttest_ind(S, SR)
    p_signif = get_pvalue_signif(p_value)

    ax = axs[idx]

    # Create boxplot
    bp = ax.boxplot([S, SR], tick_labels=['S', 'SR'], patch_artist=True)

    # Color the boxes
    colors = ['coral', 'steelblue']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Set median line color to black
    for median in bp['medians']:
        median.set_color('black')

    # Adjust y-axis with 10% padding
    all_data = np.concatenate([S, SR])
    if len(all_data) > 0:
        data_min = np.min(all_data)
        data_max = np.max(all_data)
        padding = (data_max - data_min) * 0.1
        ax.set_ylim(data_min - padding, data_max + padding)

    ax.set_title(f'{compound_name}\n{p_signif}')
    ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
out_path = os.path.join(cfg.OUTPUT_DIR, "plots", "targeted_boxplots.png")
os.makedirs(os.path.join(cfg.OUTPUT_DIR, "plots"), exist_ok=True)
plt.savefig(out_path, dpi=200, bbox_inches="tight")
#plt.show()
