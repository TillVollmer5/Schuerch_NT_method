"""
hca_dendrogram.py - Step 5b of the GCMS processing pipeline.

Produces a standalone interactive HTML dendrogram for compound class exploration:

  - Hover over any branch junction → pie-chart tooltip showing class composition
    of the subtree below that node (one pie per HCA_CLASS_ANNOTATION_COLUMNS entry)
  - Click a junction dot to add it to the comparison panel below the plot
  - Comparison panel shows selected nodes side by side with pies + legend tables
  - Click ✕ on a card (or click the dot again) to deselect
  - Colored annotation strips below the leaf axis, one row per annotation column
  - Legend entries are clickable to show / hide individual class values

All colors are consistent between pies, strips, and legends — same class always
gets the same color everywhere.

Input:  output/peak_matrix_processed.csv
        output/feature_metadata_enriched.csv  (optional; needed for strips + pies)
        output/feature_name_map.csv            (optional; FEATURE_LABEL='name')
Output: output/plots/hca_dendrogram.html

Usage:
    python hca_dendrogram.py
"""

import io
import json
import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import base64
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import plotly.graph_objects as go

import config


# --------------------------------------------------------------------------- #
#  Color palette  (shared between strips, pies, and JS comparison panel)      #
# --------------------------------------------------------------------------- #

_TAB20 = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
    "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
    "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
    "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
]
_GRAY_UNKNOWN = "#cccccc"


def _assign_colors(all_values, class_colors=None):
    """
    Map unique non-null values → hex colors.

    Pinned colors from class_colors (config.CLASS_COLORS) are applied first;
    remaining values are auto-assigned from _TAB20 in sorted alphabetical order
    so the same value always gets the same auto-assigned color.
    Unknown/NaN → _GRAY_UNKNOWN.
    """
    if class_colors is None:
        class_colors = {}
    known, seen = [], set()
    for v in all_values:
        sv = str(v).strip()
        if sv not in ("", "nan", "Unknown", "None") and sv not in seen:
            known.append(sv); seen.add(sv)
    known_sorted = sorted(known)
    auto_vals = [v for v in known_sorted if v not in class_colors]
    auto_map  = {v: _TAB20[i % len(_TAB20)] for i, v in enumerate(auto_vals)}
    cmap = {}
    for v in known_sorted:
        cmap[v] = class_colors.get(v, auto_map.get(v, _GRAY_UNKNOWN))
    cmap["Unknown"] = class_colors.get("Unknown", _GRAY_UNKNOWN)
    return cmap


# --------------------------------------------------------------------------- #
#  Linkage helpers                                                             #
# --------------------------------------------------------------------------- #

def _build_node_leaves(Z, n_leaves):
    """
    Build {node_id: [leaf_indices]} bottom-up from linkage matrix Z.
    Leaves 0..n_leaves-1; internal nodes n_leaves..2*n_leaves-2.
    """
    node_leaves = {i: [i] for i in range(n_leaves)}
    for k in range(n_leaves - 1):
        nid              = n_leaves + k
        node_leaves[nid] = node_leaves[int(Z[k, 0])] + node_leaves[int(Z[k, 1])]
    return node_leaves


# --------------------------------------------------------------------------- #
#  Pie rendering                                                               #
# --------------------------------------------------------------------------- #

def _hex_to_rgb01(hex_color):
    h = hex_color.lstrip('#')
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def _pie_b64(counts_dict, color_map, size_px=160):
    """
    Render a pie chart from {label: count} and return base64-encoded PNG.
    Transparent background so it blends into both white tooltip and white card.
    """
    labels = list(counts_dict.keys())
    values = list(counts_dict.values())
    colors = [_hex_to_rgb01(color_map.get(lbl, _GRAY_UNKNOWN)) for lbl in labels]

    inches = size_px / 100
    fig, ax = plt.subplots(figsize=(inches, inches), dpi=100)
    ax.pie(values, colors=colors, startangle=90,
           wedgeprops=dict(linewidth=0.5, edgecolor='white'))
    ax.set_aspect('equal')

    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight',
                transparent=True, dpi=100)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('ascii')


def _render_node_data(branch_fids, enriched, ann_cols, color_maps, dist):
    """
    For one internal node, compute pie images and legend tables for all
    annotation columns.  Returns a JSON-serialisable dict that is embedded
    in the HTML for the JS comparison panel.
    """
    n    = len(branch_fids)
    pies    = {}
    legends = {}   # {col: [[val, count, pct_int], ...]}

    for col in ann_cols:
        col_vals = (enriched[col]
                    .reindex(branch_fids)
                    .fillna("Unknown")
                    .astype(str)
                    .replace("nan", "Unknown"))
        counts = col_vals.value_counts().sort_values(ascending=False)

        pies[col]    = _pie_b64(counts.to_dict(), color_maps[col], size_px=160)
        legends[col] = [
            [val, int(cnt), round(100 * cnt / n)]
            for val, cnt in counts.items()
        ]

    return {"n": n, "dist": round(float(dist), 3), "pies": pies, "legends": legends}


def _tooltip_html(nd, ann_cols, color_maps):
    """
    Build the HTML string shown on hover, using pre-rendered pie images.
    One pie (+ compact legend) per annotation column, laid out side by side.
    """
    if not ann_cols or not nd.get("pies"):
        return f"<b>{nd['n']} features in subtree</b>"

    cells = []
    for col in ann_cols:
        if col not in nd["pies"]:
            continue
        b64 = nd["pies"][col]
        img = f'<img src="data:image/png;base64,{b64}" width="150" height="150">'

        rows = []
        for val, cnt, pct in nd["legends"].get(col, []):
            color  = color_maps[col].get(val, _GRAY_UNKNOWN)
            swatch = (
                f'<span style="display:inline-block;width:9px;height:9px;'
                f'background:{color};border:1px solid #aaa;'
                f'vertical-align:middle;margin-right:3px;"></span>'
            )
            rows.append(
                f'<tr><td>{swatch}</td>'
                f'<td style="font-size:9px;padding-right:5px;">{val}</td>'
                f'<td style="font-size:9px;color:#666;">{cnt}&nbsp;({pct}%)</td></tr>'
            )
        legend = (
            '<table style="border-collapse:collapse;margin-top:2px;">'
            + "".join(rows) + "</table>"
        )
        cells.append(
            f'<td style="padding:0 8px;vertical-align:top;text-align:center;">'
            f'<b style="font-size:10px;">{col}</b><br>{img}{legend}</td>'
        )

    header = f"<b>{nd['n']} features in subtree</b><br>"
    table  = '<table><tr>' + "".join(cells) + '</tr></table>'
    return header + table


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

def _build_html(plot_div, node_data_js, color_maps, ann_cols_all, ann_cols_default):
    """
    Assemble the final self-contained HTML.

    ann_cols_all     : all annotation columns with pre-computed pies
    ann_cols_default : columns pre-checked in the selector (= HCA_CLASS_ANNOTATION_COLUMNS)

    The column-selector panel lets the user toggle which columns appear in hover
    tooltips and comparison cards.  Changes call rebuildTooltips() (updates
    Plotly trace text) and rebuildPanel() (rebuilds open comparison cards).

    No str.format() is used — JS curly braces are never mis-interpreted.
    """
    node_data_json        = json.dumps(node_data_js)
    color_maps_json       = json.dumps({col: dict(cmap) for col, cmap in color_maps.items()})
    ann_cols_all_json     = json.dumps(ann_cols_all)
    ann_cols_default_json = json.dumps(ann_cols_default)

    # Checkbox panel — built with f-strings (Python values only, no JS braces)
    if ann_cols_all:
        cb_items = []
        for col in ann_cols_all:
            safe_id  = "cb_" + "".join(c if c.isalnum() else "_" for c in col)
            checked  = " checked" if col in ann_cols_default else ""
            cb_items.append(
                f'<label style="margin-right:16px;font-size:12px;cursor:pointer;">'
                f'<input type="checkbox" id="{safe_id}" onchange="onColChange()"{checked}>'
                f'&nbsp;{col}</label>'
            )
        checkbox_html = (
            '<div id="col-selector">'
            '<span style="font-weight:bold;font-size:12px;margin-right:14px;color:#444;">'
            'Pie chart columns:</span>'
            + "".join(cb_items)
            + "</div>"
        )
    else:
        checkbox_html = ""

    # JS block — plain string concatenation only, never .format()
    js = (
        "const NODE_DATA       = " + node_data_json        + ";\n"
        "const COLOR_MAPS      = " + color_maps_json       + ";\n"
        "const ANN_COLS_ALL    = " + ann_cols_all_json     + ";\n"
        "const ANN_COLS_DEFAULT= " + ann_cols_default_json + ";\n"
        "\n"
        "var selectedNodes = new Set();\n"
        "var PLOT_DIV = document.getElementById('hca-plot');\n"
        "\n"
        "function getActiveCols() {\n"
        "    return ANN_COLS_ALL.filter(function(col) {\n"
        "        var safe = '';\n"
        "        for (var i = 0; i < col.length; i++) {\n"
        "            var c = col[i];\n"
        "            safe += /[a-zA-Z0-9]/.test(c) ? c : '_';\n"
        "        }\n"
        "        var cb = document.getElementById('cb_' + safe);\n"
        "        return cb && cb.checked;\n"
        "    });\n"
        "}\n"
        "\n"
        "function buildTooltipHtml(nd, activeCols) {\n"
        "    if (!activeCols.length || !nd.pies) {\n"
        "        return '<b>' + nd.n + ' features in subtree</b>';\n"
        "    }\n"
        "    var cells = '';\n"
        "    for (var ci = 0; ci < activeCols.length; ci++) {\n"
        "        var col = activeCols[ci];\n"
        "        if (!nd.pies[col]) continue;\n"
        "        var b64 = nd.pies[col];\n"
        "        var img = '<img src=\"data:image/png;base64,' + b64 + '\" width=\"150\" height=\"150\">';\n"
        "        var rows = '';\n"
        "        var legData = nd.legends[col] || [];\n"
        "        for (var ri = 0; ri < legData.length; ri++) {\n"
        "            var val = legData[ri][0], cnt = legData[ri][1], pct = legData[ri][2];\n"
        "            var color = (COLOR_MAPS[col] && COLOR_MAPS[col][val]) || '#cccccc';\n"
        "            var swatch = '<span style=\"display:inline-block;width:9px;height:9px;'\n"
        "                + 'background:' + color + ';border:1px solid #aaa;'\n"
        "                + 'vertical-align:middle;margin-right:3px;\"></span>';\n"
        "            rows += '<tr><td>' + swatch + '</td>'\n"
        "                + '<td style=\"font-size:9px;padding-right:5px;\">' + val + '</td>'\n"
        "                + '<td style=\"font-size:9px;color:#666;\">' + cnt + '&nbsp;(' + pct + '%)</td></tr>';\n"
        "        }\n"
        "        var legend = '<table style=\"border-collapse:collapse;margin-top:2px;\">' + rows + '</table>';\n"
        "        cells += '<td style=\"padding:0 8px;vertical-align:top;text-align:center;\">'\n"
        "            + '<b style=\"font-size:10px;\">' + col + '</b><br>' + img + legend + '</td>';\n"
        "    }\n"
        "    var header = '<b>' + nd.n + ' features in subtree</b><br>';\n"
        "    var table = '<table><tr>' + cells + '</tr></table>';\n"
        "    return header + table;\n"
        "}\n"
        "\n"
        "function rebuildTooltips() {\n"
        "    var activeCols = getActiveCols();\n"
        "    var idx = getNodeTraceIdx();\n"
        "    if (idx < 0) return;\n"
        "    var n = PLOT_DIV.data[idx].x.length;\n"
        "    var newTexts = [];\n"
        "    for (var k = 0; k < n; k++) {\n"
        "        newTexts.push(buildTooltipHtml(NODE_DATA[k], activeCols));\n"
        "    }\n"
        "    Plotly.restyle('hca-plot', {'text': [newTexts]}, [idx]);\n"
        "}\n"
        "\n"
        "function onColChange() {\n"
        "    rebuildTooltips();\n"
        "    rebuildPanel();\n"
        "}\n"
        "\n"
        "PLOT_DIV.on('plotly_click', function(data) {\n"
        "    if (!data.points || data.points.length === 0) return;\n"
        "    var pt = data.points[0];\n"
        "    if (pt.data.name !== '__node_markers__') return;\n"
        "    var k = pt.pointIndex;\n"
        "    if (selectedNodes.has(k)) { selectedNodes.delete(k); }\n"
        "    else { selectedNodes.add(k); }\n"
        "    updateMarkers();\n"
        "    rebuildPanel();\n"
        "});\n"
        "\n"
        "function getNodeTraceIdx() {\n"
        "    return PLOT_DIV.data.findIndex(function(t) { return t.name === '__node_markers__'; });\n"
        "}\n"
        "\n"
        "function updateMarkers() {\n"
        "    var idx = getNodeTraceIdx();\n"
        "    if (idx < 0) return;\n"
        "    var n = PLOT_DIV.data[idx].x.length;\n"
        "    var colors = [], sizes = [];\n"
        "    for (var k = 0; k < n; k++) {\n"
        "        colors.push(selectedNodes.has(k) ? '#e74c3c' : 'rgba(70,70,70,0.30)');\n"
        "        sizes.push(selectedNodes.has(k) ? 13 : 9);\n"
        "    }\n"
        "    Plotly.restyle('hca-plot', {'marker.color': [colors], 'marker.size': [sizes]}, [idx]);\n"
        "}\n"
        "\n"
        "function rebuildPanel() {\n"
        "    var panel = document.getElementById('comparison-panel');\n"
        "    var cards = document.getElementById('comparison-cards');\n"
        "    if (selectedNodes.size === 0) { panel.style.display = 'none'; return; }\n"
        "    panel.style.display = 'block';\n"
        "    var html = '';\n"
        "    var sortedNodes = Array.from(selectedNodes).sort(function(a, b) { return a - b; });\n"
        "    for (var i = 0; i < sortedNodes.length; i++) {\n"
        "        html += buildCard(sortedNodes[i], NODE_DATA[sortedNodes[i]]);\n"
        "    }\n"
        "    cards.innerHTML = html;\n"
        "}\n"
        "\n"
        "function buildCard(k, nd) {\n"
        "    var activeCols = getActiveCols();\n"
        "    var piesHtml = '';\n"
        "    for (var ci = 0; ci < activeCols.length; ci++) {\n"
        "        var col = activeCols[ci];\n"
        "        if (!nd.pies || !nd.pies[col]) continue;\n"
        "        var legendRows = '';\n"
        "        var legData = nd.legends[col] || [];\n"
        "        for (var ri = 0; ri < legData.length; ri++) {\n"
        "            var val = legData[ri][0], cnt = legData[ri][1], pct = legData[ri][2];\n"
        "            var color = (COLOR_MAPS[col] && COLOR_MAPS[col][val]) || '#cccccc';\n"
        "            legendRows += '<tr>'\n"
        "                + '<td style=\"padding:1px 2px 1px 0;\">'\n"
        "                + '<span style=\"display:inline-block;width:9px;height:9px;'\n"
        "                + 'background:' + color + ';border:1px solid #aaa;'\n"
        "                + 'vertical-align:middle;\"></span></td>'\n"
        "                + '<td style=\"font-size:9px;padding-right:6px;vertical-align:middle;\">'\n"
        "                + val + '</td>'\n"
        "                + '<td style=\"font-size:9px;color:#666;vertical-align:middle;\">'\n"
        "                + cnt + ' (' + pct + '%)</td></tr>';\n"
        "        }\n"
        "        var legendTable = '<table style=\"border-collapse:collapse;margin-top:3px;\">'\n"
        "            + legendRows + '</table>';\n"
        "        piesHtml += '<div style=\"text-align:center;margin:0 12px 0 0;\">'\n"
        "            + '<div style=\"font-weight:bold;font-size:11px;margin-bottom:3px;\">'\n"
        "            + col + '</div>'\n"
        "            + '<img src=\"data:image/png;base64,' + nd.pies[col] + '\"'\n"
        "            + ' width=\"160\" height=\"160\"'\n"
        "            + ' style=\"display:block;margin:0 auto;\"/>'\n"
        "            + legendTable + '</div>';\n"
        "    }\n"
        "    if (!piesHtml) {\n"
        "        piesHtml = '<span style=\"font-size:11px;color:#888;\">'\n"
        "            + 'No class annotation available.</span>';\n"
        "    }\n"
        "    return '<div style=\"border:1px solid #ddd;border-radius:6px;padding:12px;'\n"
        "        + 'background:#fff;position:relative;min-width:180px;\">'\n"
        "        + '<div style=\"font-weight:bold;font-size:12px;margin-bottom:8px;'\n"
        "        + 'padding-right:24px;\">Node ' + (k+1)\n"
        "        + '<span style=\"font-weight:normal;color:#555;font-size:11px;\">'\n"
        "        + ' \u2014 ' + nd.n + ' features&nbsp; dist: ' + nd.dist + '</span></div>'\n"
        "        + '<button onclick=\"removeNode(' + k + ')\"'\n"
        "        + ' style=\"position:absolute;top:8px;right:8px;border:none;background:none;'\n"
        "        + 'cursor:pointer;font-size:16px;color:#bbb;line-height:1;\">\u2715</button>'\n"
        "        + '<div style=\"display:flex;flex-wrap:wrap;\">'\n"
        "        + piesHtml + '</div></div>';\n"
        "}\n"
        "\n"
        "function removeNode(k) {\n"
        "    selectedNodes.delete(k);\n"
        "    updateMarkers();\n"
        "    rebuildPanel();\n"
        "}\n"
    )

    parts = [
        "<!DOCTYPE html>\n<html>\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<title>HCA Class Dendrogram</title>\n"
        "<style>\n"
        "  body  { font-family: Arial, sans-serif; margin: 0; padding: 0; background: #fff; }\n"
        "  #col-selector {\n"
        "    display: flex; align-items: center; flex-wrap: wrap;\n"
        "    padding: 10px 20px; background: #f5f5f5;\n"
        "    border-bottom: 1px solid #ddd;\n"
        "  }\n"
        "  #comparison-panel {\n"
        "    display: none; padding: 14px 20px 20px;\n"
        "    border-top: 2px solid #e0e0e0; background: #f8f8f8;\n"
        "  }\n"
        "  #comparison-panel h3 { margin: 0 0 4px; font-size: 14px; color: #333; }\n"
        "  .hint { color: #999; font-size: 11px; margin-bottom: 12px; }\n"
        "  #comparison-cards {\n"
        "    display: flex; flex-wrap: wrap; gap: 16px;\n"
        "    max-height: 680px; overflow-y: auto;\n"
        "  }\n"
        "</style>\n"
        "</head>\n<body>\n",

        checkbox_html,

        plot_div,

        "\n<div id=\"comparison-panel\">\n"
        "  <h3>Node comparison</h3>\n"
        "  <div class=\"hint\">\n"
        "    Click junction dots in the dendrogram to add nodes here. "
        "Click \u2715 on a card (or click the dot again) to remove it.\n"
        "  </div>\n"
        "  <div id=\"comparison-cards\"></div>\n"
        "</div>\n"
        "<script>\n",

        js,

        "</script>\n</body>\n</html>\n",
    ]
    return "".join(parts)


def run(cfg=config):
    if not getattr(cfg, "RUN_HCA_DENDROGRAM", True):
        return

    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("-- Step 5b: Interactive class dendrogram -------------------------")

    matrix_path   = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed.csv")
    enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
    name_map_path = os.path.join(cfg.OUTPUT_DIR, "feature_name_map.csv")

    if not os.path.exists(matrix_path):
        print("  [skip] peak_matrix_processed.csv not found; run normalization.py first")
        return

    matrix      = pd.read_csv(matrix_path, index_col="sample")
    feature_ids = list(matrix.columns)
    n_features  = len(feature_ids)
    print(f"  features : {n_features}")

    # --- enriched metadata ---------------------------------------------------
    ann_cols_default = list(getattr(cfg, "HCA_CLASS_ANNOTATION_COLUMNS", []))
    class_colors     = getattr(cfg, "CLASS_COLORS", {})

    enriched = None
    if os.path.exists(enriched_path):
        enriched = pd.read_csv(enriched_path, index_col="feature_id")
        ann_cols_default = [c for c in ann_cols_default if c in enriched.columns]
    else:
        ann_cols_default = []
        print("  [info] feature_metadata_enriched.csv not found — plain dendrogram")

    # ann_cols = columns with pre-computed pies (= configured columns)
    # Same as default here; users expand by editing HCA_CLASS_ANNOTATION_COLUMNS.
    ann_cols = ann_cols_default

    # --- color maps (built from full dataset so colors are globally consistent)
    color_maps = {}
    for col in ann_cols:
        all_vals = (enriched[col].fillna("Unknown").astype(str)
                    .replace("nan", "Unknown").tolist())
        color_maps[col] = _assign_colors(all_vals, class_colors)

    # --- feature display labels ----------------------------------------------
    label_map = {}
    if getattr(cfg, "FEATURE_LABEL", "id") == "name" and os.path.exists(name_map_path):
        nm = pd.read_csv(name_map_path, index_col="feature_id")
        if "compound_name" in nm.columns:
            label_map = {
                fid: nm.at[fid, "compound_name"]
                for fid in feature_ids
                if fid in nm.index
                and pd.notna(nm.at[fid, "compound_name"])
                and str(nm.at[fid, "compound_name"]).strip()
            }

    label_col = getattr(cfg, "CLASS_LABEL_COLUMN", "")

    def _display(fid):
        base = label_map.get(fid, fid)
        if label_col and enriched is not None and label_col in enriched.columns:
            if fid in enriched.index:
                v = enriched.at[fid, label_col]
                if pd.notna(v) and str(v).strip() not in ("", "nan"):
                    return f"{base} [{v}]"
        return base

    # --- linkage -------------------------------------------------------------
    linkage_method = getattr(cfg, "HCA_LINKAGE", "ward")
    metric         = getattr(cfg, "HCA_METRIC",  "euclidean")

    data_arr = matrix.values.T
    dist_mat = ssd.pdist(data_arr, metric=metric)
    Z        = sch.linkage(dist_mat, method=linkage_method)
    max_dist = float(Z[-1, 2])

    dend_data        = sch.dendrogram(Z, no_plot=True)
    ordered_leaf_idx = dend_data['leaves']
    ordered_fids     = [feature_ids[i] for i in ordered_leaf_idx]
    n_int            = n_features - 1

    leaf_xpos = {orig_idx: 5 + 10 * rank
                 for rank, orig_idx in enumerate(ordered_leaf_idx)}

    node_leaves = _build_node_leaves(Z, n_features)

    node_xpos = dict(leaf_xpos)
    for k in range(n_int):
        nid            = n_features + k
        node_xpos[nid] = (node_xpos[int(Z[k, 0])] + node_xpos[int(Z[k, 1])]) / 2

    # --- build node data (pie images + legend tables) for every internal node -
    print(f"  rendering node pies : ", end="", flush=True)
    node_data_js = {}   # k (0-indexed merge order) → JSON-serialisable dict

    hov_x, hov_y, hov_html = [], [], []

    for k in range(n_int):
        nid         = n_features + k
        branch_idxs = node_leaves[nid]
        branch_fids = [feature_ids[i] for i in branch_idxs]

        if enriched is not None and ann_cols:
            nd = _render_node_data(
                branch_fids, enriched, ann_cols, color_maps, float(Z[k, 2])
            )
        else:
            nd = {"n": len(branch_fids), "dist": round(float(Z[k, 2]), 3),
                  "pies": {}, "legends": {}}

        node_data_js[k] = nd
        hov_x.append(node_xpos[nid])
        hov_y.append(float(Z[k, 2]))
        hov_html.append(_tooltip_html(nd, ann_cols, color_maps))

        if (k + 1) % 25 == 0 or k == n_int - 1:
            print(f"{k+1}", end=" ", flush=True)

    print()

    # --- traces: dendrogram lines --------------------------------------------
    line_x, line_y = [], []
    for xs, ys in zip(dend_data['icoord'], dend_data['dcoord']):
        line_x.extend(list(xs) + [None])
        line_y.extend(list(ys) + [None])

    traces = [
        go.Scatter(
            x=line_x, y=line_y,
            mode='lines',
            line=dict(color='#555555', width=1.2),
            hoverinfo='skip',
            showlegend=False,
        )
    ]

    # --- trace: clickable node markers (hover tooltip + click to select) -----
    traces.append(
        go.Scatter(
            x=hov_x, y=hov_y,
            mode='markers',
            marker=dict(
                size=9,
                color='rgba(70,70,70,0.30)',
                symbol='circle',
                line=dict(width=1, color='rgba(70,70,70,0.55)'),
            ),
            hovertemplate="%{text}<extra></extra>",
            text=hov_html,
            showlegend=False,
            name='__node_markers__',
            hoverlabel=dict(
                bgcolor='white',
                bordercolor='#cccccc',
                font=dict(size=10),
                namelength=0,
            ),
        )
    )

    # --- traces: annotation strips below y=0 ---------------------------------
    strip_dy         = max_dist * 0.09
    n_strips         = 0
    group_title_done = set()
    marker_px        = max(5, min(14, int(700 / n_features)))

    if enriched is not None and ann_cols:
        n_strips = len(ann_cols)
        for ci, col in enumerate(ann_cols):
            strip_y = -(ci + 0.5) * strip_dy
            cmap    = color_maps[col]
            vals    = (enriched[col]
                       .reindex(ordered_fids)
                       .fillna("Unknown")
                       .astype(str)
                       .replace("nan", "Unknown"))

            for val, color in cmap.items():
                positions = [5 + 10 * pos for pos, v in enumerate(vals) if v == val]
                fids_here = [ordered_fids[pos] for pos, v in enumerate(vals) if v == val]
                if not positions:
                    continue
                first_in_group = col not in group_title_done
                group_title_done.add(col)
                traces.append(go.Scatter(
                    x=positions, y=[strip_y] * len(positions),
                    mode='markers',
                    marker=dict(symbol='square', size=marker_px,
                                color=color, line=dict(width=0)),
                    name=val,
                    legendgroup=col,
                    legendgrouptitle_text=col if first_in_group else None,
                    showlegend=True,
                    hovertemplate=(
                        f"<b>{col}:</b> {val}<br>"
                        "Feature: %{customdata}<extra></extra>"
                    ),
                    customdata=fids_here,
                ))

    # --- layout --------------------------------------------------------------
    y_min   = -(n_strips + 0.3) * strip_dy if n_strips else -max_dist * 0.04
    fig_w   = max(900,  min(n_features * 15 + 260, 4200))
    fig_h   = max(520,  560 + n_strips * 45)

    tick_vals = [5 + 10 * i for i in range(n_features)]
    tick_text = [_display(fid) for fid in ordered_fids]

    annotations = []
    for ci, col in enumerate(ann_cols):
        annotations.append(dict(
            x=-18, y=-(ci + 0.5) * strip_dy,
            text=f"<i>{col}</i>",
            xanchor='right', xref='x', yref='y',
            showarrow=False, font=dict(size=8, color='#444'),
        ))

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(
            text=(
                f"Interactive compound dendrogram — "
                f"{cfg.HCA_LINKAGE} / {cfg.HCA_METRIC} | {n_features} features<br>"
                "<sup>Hover junction dots for class pies  ·  "
                "Click a dot to add it to the comparison panel below</sup>"
            ),
            x=0.5, font=dict(size=13),
        ),
        xaxis=dict(
            tickmode='array', tickvals=tick_vals, ticktext=tick_text,
            tickangle=-90,
            tickfont=dict(size=max(6, min(9, int(600 / n_features)))),
            range=[-35, 10 * n_features + 25],
            showgrid=False, zeroline=False,
        ),
        yaxis=dict(
            title='Linkage distance',
            range=[y_min, max_dist * 1.06],
            showgrid=True, gridcolor='#f0f0f0',
            zeroline=True, zerolinecolor='#aaaaaa', zerolinewidth=1,
        ),
        plot_bgcolor='white', paper_bgcolor='white',
        width=fig_w, height=fig_h,
        hovermode='closest',
        legend=dict(
            tracegroupgap=10, font=dict(size=9),
            title=dict(text='<b>Compound classes</b>', font=dict(size=10)),
            bgcolor='rgba(255,255,255,0.92)',
            bordercolor='#dddddd', borderwidth=1,
        ),
        margin=dict(l=90, r=230, t=85, b=160),
        annotations=annotations,
    )

    # --- write self-contained HTML -------------------------------------------
    plot_div = fig.to_html(
        full_html=False,
        include_plotlyjs='cdn',
        div_id='hca-plot',
        config={'responsive': False},
    )

    html = _build_html(plot_div, node_data_js, color_maps, ann_cols, ann_cols_default)

    out_path = os.path.join(plots_dir, "hca_dendrogram.html")
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"  -> {out_path}")
    if n_strips:
        print(f"     annotation strips : {ann_cols}")
    print("     Hover dots → pie tooltip  |  Click dot → add to comparison panel")

    return fig


if __name__ == "__main__":
    run()
