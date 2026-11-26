#!/usr/bin/env python3

"""
This script can be used to generate a heatmap with the cg/wgMLST alleles distribution across a species tree

By Vitor Borges
@INSA
"""

import sys
import argparse
import datetime as datetime
import pandas as pd
import numpy as np
from Bio import Phylo
import colorsys
from collections import defaultdict
import plotly.graph_objects as go

pd.set_option('future.no_silent_downcasting', True)

version = "2.0.0"
last_updated = "2026-11-26"

SAMPLE_INTEREST_COLOR = "#800080"

def print_log(message, log):
	""" print messages in the terminal and in the log file """

	print(message)
	print(message, file = log)

def lighten_color(hex_color, factor=0.5):
    """
    get fainted colors
    """
    
    hex_color = hex_color.strip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    h, l, s = colorsys.rgb_to_hls(r / 255, g / 255, b / 255)
    l = min(1, l + (1 - l) * factor)
    r, g, b = colorsys.hls_to_rgb(h, l, s)

    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"

def generate_auto_colors(categories):
    """
    Automatically generate colours for metadata categories.
    """

    cats = list(categories)
    n = max(1, len(cats))
    color_map_cat = {}

    for i, cat in enumerate(cats):
        h = i / n
        s = 0.6
        v = 0.8
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        color_map_cat[cat] = (
            f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
        )

    return color_map_cat

def read_tree_order(tree_file):
    """
    get a list with the samples of a tree
    """

    tree = Phylo.read(tree_file, "newick")

    return [term.name for term in tree.get_terminals()]

def read_sample_of_interest(sample_file):
    """
    read sample IDs from the first column of a TSV file.
    """

    df = pd.read_csv(sample_file, sep="\t", header=None, dtype=str)
    samples = set(df.iloc[:, 0].dropna().astype(str).str.strip())

    return set(samples)

def assign_colors(allelic_matrix, metadata, meta_col, color_map, mode="frequency", ncat=None, ncat_dominant=None, sample_of_interest=None, sample_interest_color=SAMPLE_INTEREST_COLOR):
    """
    assign colors and classify alleles.

    IMPORTANT:
    Samples of interest are completely excluded from category-based
    analyses. Their metadata category is ignored.

    They are considered separately only for identifying alleles that
    occur exclusively in the sample(s) of interest.
    """

    sample_of_interest = set(sample_of_interest or [])

    allele_to_color = {}

    exclusives = defaultdict(list)
    shared = defaultdict(list)

    # Alleles exclusively found in samples of interest
    sample_interest_exclusives = []

    allele_to_categories = defaultdict(set)

    conserved_loci = set()

    # Only samples that are not samples of interest are used for category calculations
    regular_metadata = metadata.loc[~metadata.index.isin(sample_of_interest)]

    all_cats = regular_metadata[meta_col].dropna().unique()
    total_per_cat = regular_metadata[meta_col].value_counts()

    for locus in allelic_matrix.columns:
        locus_data = allelic_matrix[locus]
        locus_data = locus_data.replace("0", np.nan).infer_objects(copy=False)
        unique_alleles = locus_data.dropna().unique()

        regular_locus_data = locus_data.loc[~locus_data.index.isin(sample_of_interest)].dropna()

        if len(unique_alleles) == 0:
            continue
        elif len(unique_alleles) == 1 and len(regular_locus_data) > 0:
            conserved_loci.add(locus)
            only_allele = unique_alleles[0]
            allele_to_color[(locus, only_allele)] = "#000000"
            allele_to_categories[(locus, only_allele)] = set(all_cats)
        else:
            for allele in unique_alleles:
                idx_all = locus_data[locus_data == allele].index
                soi_idx = idx_all[idx_all.isin(sample_of_interest)]
                regular_idx = idx_all[~idx_all.isin(sample_of_interest)]

                categories = metadata.loc[regular_idx, meta_col].dropna()
                counts = categories.value_counts()
                category_set = set(categories.unique())
                allele_to_categories[(locus, allele)] = category_set

                # CASE 1: Allele occurs in one or more categories

                if len(category_set) > 0:
                    if mode == "frequency":
                        freqs = counts / total_per_cat[counts.index]
                        dominant_cat = freqs.idxmax()
                    else:
                        dominant_cat = counts.idxmax()

                    if len(category_set) == 1:
                        base_color = color_map[dominant_cat]
                        allele_to_color[(locus, allele)] = lighten_color(base_color,factor=0.5)
                        exclusives[dominant_cat].append((locus, allele))
                    else:
                        base_color = color_map[dominant_cat]
                        allele_to_color[(locus, allele)] = base_color
                        shared[(locus, allele)] = sorted(category_set)

                    if ncat is not None and len(category_set) >= ncat:
                        allele_to_color[(locus, allele)] = "#000000"

                    if ncat_dominant is not None:
                        dominant_count = 0
                        for cat in all_cats:
                            cat_idx = regular_metadata[regular_metadata[meta_col] == cat].index

                            if len(cat_idx) == 0:
                                continue
                            cat_series = locus_data.loc[locus_data.index.intersection(cat_idx)]
                            if len(cat_series.dropna()) == 0:
                                continue
                            allele_counts = cat_series.value_counts()
                            if allele_counts.idxmax() == allele:
                                dominant_count += 1
                        if dominant_count >= ncat_dominant:
                            allele_to_color[(locus, allele)] = "#000000"

                # CASE 2: Allele occurs ONLY in samples of interest
                else:
                    if len(soi_idx) > 1:
                        soi_type = "exclusive (group sample_of_interest)"
                    elif len(soi_idx) == 1:
                        soi_type = "exclusive (this sample_of_interest)"
                    else:
                        continue

                    allele_to_color[(locus, allele)] = sample_interest_color
                    sample_interest_exclusives.append((locus, allele, soi_type))

    return allele_to_color, exclusives, shared, sample_interest_exclusives, allele_to_categories, conserved_loci

def build_heatmap(log, allelic_matrix, allele_to_color, allele_to_categories, sample_order, color_map, meta_col, output_prefix, sample_of_interest=None, sample_interest_exclusives=None, tree_file=None, tree_labels=None, tree_scale=None, sample_interest_color=SAMPLE_INTEREST_COLOR):
    """
    generate the heatmap
    """

    sample_of_interest = set(sample_of_interest or [])
    sample_interest_exclusives = sample_interest_exclusives or []

    sample_interest_exclusive_map = {(locus, allele): soi_type for locus, allele, soi_type in sample_interest_exclusives}
    loci = allelic_matrix.columns.tolist()
    samples = [s for s in sample_order if s in allelic_matrix.index]

    unique_colors = list({v for v in allele_to_color.values()})
    if "#ffffff" not in unique_colors:
        unique_colors.append("#ffffff")
    if "#000000" not in unique_colors:
        unique_colors.append("#000000")
    if sample_of_interest and sample_interest_color not in unique_colors:
        unique_colors.append(sample_interest_color)

    color_to_id = {color: i for i, color in enumerate(unique_colors)}

    z = []
    hover_text = []
    for sample in samples:
        row = []
        h_row = []
        for locus in loci:
            allele = allelic_matrix.loc[sample, locus]
            if allele == "0" or pd.isna(allele):
                row.append(color_to_id["#ffffff"])
                h_row.append(f"Sample: {sample}<br>" f"{locus}<br>" f"Allele: 0<br>" f"Category: -")
            else:
                color = allele_to_color.get((locus, allele), "#ffffff")
                row.append(color_to_id[color])
                cats = allele_to_categories.get((locus, allele), set())

                soi_type = sample_interest_exclusive_map.get((locus, allele))
                if soi_type is not None:
                    cat_str = soi_type
                elif cats:
                    cat_str = ", ".join(sorted(list(cats)))
                else:
                    cat_str = "-"

                h_row.append(f"Sample: {sample}<br>" f"{locus}<br>" f"Allele: {allele}<br>" f"Category: {cat_str}")
        z.append(row)
        hover_text.append(h_row)

    fig = go.Figure()
    sample_y = {sample: i for i, sample in enumerate(samples)}
    heatmap_y = [sample_y[sample] for sample in samples]
    fig.add_trace(go.Heatmap(z=z,x=[f"Locus {i + 1}" for i in range(len(loci))],y=heatmap_y,hoverinfo="text",text=hover_text,colorscale=[[i / (len(unique_colors) - 1),c] for i, c in enumerate(unique_colors)], showscale=False, xaxis="x2", yaxis="y"))

    if tree_file is not None and samples:
        tree = Phylo.read(tree_file, "newick")
        depths = tree.depths()

        if not depths or max(depths.values()) == 0:
            depths = tree.depths(unit_branch_lengths=True)

        tree_x = []
        tree_y = []

        def get_clade_y(clade):
            """Return y coordinate for a clade using only samples
            that are present in the heatmap."""

            if clade.is_terminal():
                if clade.name in sample_y:
                    return float(sample_y[clade.name])
                return None

            child_y = []
            for child in clade.clades:
                y = get_clade_y(child)
                if y is not None:
                    child_y.append(y)
            if not child_y:
                return None

            y_parent = sum(child_y) / len(child_y)
            x_parent = depths.get(clade, 0.0)
            tree_x.extend([x_parent,x_parent,None])
            tree_y.extend([min(child_y),max(child_y),None])

            for child in clade.clades:
                y_child = get_clade_y_for_draw(child)
                if y_child is None:
                    continue
                x_child = depths.get(child, x_parent)
                tree_x.extend([x_parent,x_child,None])
                tree_y.extend([y_child,y_child,None])

            return y_parent

        def get_clade_y_for_draw(clade):
            if clade.is_terminal():
                if clade.name in sample_y:
                    return float(sample_y[clade.name])
                return None

            child_y = []
            for child in clade.clades:
                y = get_clade_y_for_draw(child)
                if y is not None:
                    child_y.append(y)

            if not child_y:
                return None

            return sum(child_y) / len(child_y)

        tree_x = []
        tree_y = []

        def draw_clade(clade):
            if clade.is_terminal():
                if clade.name not in sample_y:
                    return None
                return float(sample_y[clade.name])

            child_positions = []

            for child in clade.clades:
                y_child = draw_clade(child)
                if y_child is not None:
                    child_positions.append((child,y_child))

            if not child_positions:
                return None

            x_parent = depths.get(clade, 0.0)
            y_parent = sum(y for _, y in child_positions) / len(child_positions)
            tree_x.extend([x_parent,x_parent,None])
            tree_y.extend([min(y for _, y in child_positions),max(y for _, y in child_positions),None])

            for child, y_child in child_positions:
                x_child = depths.get(child,x_parent)
                tree_x.extend([x_parent,x_child,None])
                tree_y.extend([y_child,y_child,None])
            
            return y_parent

        draw_clade(tree.root)
        fig.add_trace(go.Scatter(x=tree_x,y=tree_y,mode="lines",line=dict(color="#444444",width=1.5),hoverinfo="skip",showlegend=False,xaxis="x",yaxis="y"))

        if tree_labels:
            tip_samples = []
            tip_x = []
            tip_y = []
            for terminal in tree.get_terminals():
                if terminal.name in sample_y:
                    tip_samples.append(terminal.name)
                    tip_x.append(depths.get(terminal, 0.0))
                    tip_y.append(sample_y[terminal.name])
            fig.add_trace(go.Scatter(x=tip_x,y=tip_y,mode="text",text=tip_samples,textposition="middle right",textfont=dict(size=tree_labels),cliponaxis=False,hoverinfo="skip",showlegend=False,xaxis="x",yaxis="y"))

    tree_max_depth = 0.0
    if tree_file is not None and samples:
        try:
            tree_max_depth = max(depths.get(terminal, 0.0) for terminal in tree.get_terminals() if terminal.name in sample_y)
        except (ValueError, TypeError):
            tree_max_depth = 0.0

    if sample_of_interest and samples:
        soi_y = [sample_y[sample] for sample in samples if sample in sample_of_interest]
        if soi_y:
            fig.add_trace(go.Scatter(x=[0.5] * len(soi_y),y=soi_y,mode="markers",marker=dict(symbol="square",size=7,color="#000000"),hoverinfo="text",text=[f"Sample of interest: {sample}" for sample in samples if sample in sample_of_interest],showlegend=False,xaxis="x3",yaxis="y"))

    fig.update_layout(
        xaxis=dict(
            domain=[0.0, 0.20],
            showticklabels=(tree_scale is not None),
            showgrid=False,
            zeroline=False,
            ticks=("outside" if tree_scale is not None else ""),
            tickfont=dict(size=tree_scale) if tree_scale is not None else None,
            showline=(tree_scale is not None),
            linewidth=1,
            mirror=False,
            tickvals=(
                [
                    0.0,
                    tree_max_depth / 4.0,
                    tree_max_depth / 2.0,
                    3.0 * tree_max_depth / 4.0,
                    tree_max_depth
               ]
                if tree_scale is not None and tree_max_depth > 0
                else None
           ),
            ticktext=(
                [
                    f"{tree_max_depth:g}",
                    f"{tree_max_depth * 0.75:g}",
                    f"{tree_max_depth * 0.5:g}",
                    f"{tree_max_depth * 0.25:g}",
                    "0"
               ]
                if tree_scale is not None and tree_max_depth > 0
                else None
           )
       ),
        xaxis3=dict(
            domain=[0.20, 0.21],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            fixedrange=True
       ),
        xaxis2=dict(
            domain=[0.21, 1.0],
            showticklabels=False,
            showgrid=False,
            zeroline=False
       ),
        yaxis=dict(
            showticklabels=False,
            autorange="reversed",
            showgrid=False,
            zeroline=False,
            fixedrange=False
       ),
        margin=dict(
            t=50,
            l=50,
            r=250
       ),
        title="Allelic Category Heatmap"
   )

    annotations = []

    legend_y = 1
    step = 0.05
    for cat, base_color in color_map.items():
        faint_color = lighten_color(base_color,0.5)
        annotations.append(
            dict(x=1.02,y=legend_y,xref="paper",yref="paper",xanchor="left",text=f"{cat} (faint=exclusive)",showarrow=False,font=dict(color=faint_color)))
        legend_y -= step
        annotations.append(dict(x=1.02,y=legend_y,xref="paper",yref="paper",xanchor="left",text=f"{cat} (dark=shared)",showarrow=False,font=dict(color=base_color)))
        legend_y -= step

    if sample_of_interest:
        annotations.append(dict(x=1.02,y=legend_y,xref="paper",yref="paper",xanchor="left",text="sample_of_interest (exclusive)",showarrow=False,font=dict(color=sample_interest_color)))

    fig.update_layout(annotations=annotations)
    html_out = f"{output_prefix}.html"

    extra_buttons = f"""
<script>
function download(format) {{
    Plotly.downloadImage(
        document.getElementsByClassName('js-plotly-plot')[0],
        {{
            format: format,
            filename: '{output_prefix}'
       }}
   )
}}
</script>

<div style="margin-bottom:10px;">
  <button onclick="download('png')">
    Download PNG
  </button>
  <button onclick="download('svg')">
    Download SVG
  </button>
</div>
"""

    with open(html_out, "w") as f:
        f.write(extra_buttons)
        f.write(fig.to_html(include_plotlyjs="cdn",full_html=False))

    print_log(f"\t[OK] Heatmap saved to {html_out}", log)

def save_summary(log, exclusives,shared,sample_interest_exclusives,allele_to_categories,output_prefix):
    """
    generate a summary report
    """
    rows = []

    for cat, alleles in exclusives.items():
        for locus, allele in alleles:
            rows.append([cat,"exclusive",locus,allele])

    for (locus, allele), cats in shared.items():
        cat_str = ",".join(sorted(cats))
        rows.append([cat_str,"shared",locus,allele])

    for locus, allele, soi_type in sample_interest_exclusives:
        rows.append(["sample_of_interest",soi_type,locus,allele])

    df = pd.DataFrame(rows,columns=["Category","Type","Locus","Allele"])
    out1 = (f"{output_prefix}_allele_summary.tsv")

    df.to_csv(out1, sep="\t", index=False)

    print_log(f"\t[OK] Allele summary saved to {out1}", log)

    # Summary counts
    if not df.empty:
        summary_counts = (df.groupby(["Category", "Type"]).agg(n_loci=("Locus", "nunique"),n_alleles=("Allele", "count")).reset_index())

    else:
        summary_counts = pd.DataFrame(columns=["Category","Type","n_loci","n_alleles"])

    out2 = (f"{output_prefix}_allele_summary_counts.tsv")

    summary_counts.to_csv(out2,sep="\t",index=False)

    print_log(f"\t[OK] Allele summary counts saved to {out2}", log)

def save_auto_colors_tsv(log, color_map, output_prefix, meta_col, sample_interest_color=SAMPLE_INTEREST_COLOR):

    rows = []

    for cat, color in color_map.items():
        rows.append([meta_col, cat, color])

    # Add the special sample_of_interest colour
    rows.append([meta_col, "sample_of_interest", sample_interest_color])

    df = pd.DataFrame(rows, columns=["column", "category", "color"])

    out = f"{output_prefix}_colors.tsv"
    df.to_csv(out, sep="\t", index=False)

    print_log(f"\n[OK] Auto-generated colors saved to {out}", log)

def build_detailed_matrix(log, allelic_matrix, sample_order, allele_to_categories, conserved_loci, output_prefix, sample_interest_exclusives=None):
    
    sample_interest_exclusives = (sample_interest_exclusives or [])
    sample_interest_exclusive_map = {(locus, allele): soi_type for locus, allele, soi_type in sample_interest_exclusives}
    loci = allelic_matrix.columns.tolist()
    samples = [s for s in sample_order if s in allelic_matrix.index]
    data = []

    for sample in samples:
        row = []
        for locus in loci:
            allele = allelic_matrix.loc[sample, locus]
            if allele == "0" or pd.isna(allele):
                row.append("0 (missing; -)")
            else:
                soi_type = (sample_interest_exclusive_map.get((locus, allele)))
                if soi_type is not None:
                    annot = (f"{allele} ({soi_type})")
                else:
                    cats = allele_to_categories.get((locus, allele), set())
                    if locus in conserved_loci:
                        annot = (f"{allele} " f"(conserved; all)")
                    else:
                        if len(cats) <= 1:
                            cat_str = ("-" if len(cats) == 0 else ", ".join(sorted(list(cats))))
                            annot = (f"{allele} " f"(exclusive; {cat_str})")
                        else:
                            cat_str = ", ".join(sorted(list(cats)))
                            annot = (f"{allele} " f"(shared; {cat_str})")
                row.append(annot)
        data.append(row)

    df = pd.DataFrame(data, index=samples, columns=loci)
    out = (f"{output_prefix}_allelic_detailed_matrix.tsv")
    df.to_csv(out, sep="\t")

    print_log(f"\t[OK] Detailed matrix saved to {out}", log)

def build_dominance_table(log, allelic_matrix, metadata, meta_col, mode, output_prefix, sample_of_interest=None):
    """
    Build the allele dominance table.

    Samples of interest are completely excluded from:
      - category counts
      - category frequencies
      - dominant category
    """

    sample_of_interest = set(sample_of_interest or [])
    regular_metadata = metadata.loc[~metadata.index.isin(sample_of_interest)]
    all_cats = list(regular_metadata[meta_col].dropna() .unique())
    total_per_cat = (regular_metadata[meta_col].value_counts())

    rows = []

    for locus in allelic_matrix.columns:
        locus_series = (allelic_matrix[locus] .replace("0", np.nan))
        present = (locus_series .dropna())
        unique_alleles = present.unique()
        for allele in unique_alleles:
            idx_all = present[present == allele].index
            soi_idx = idx_all[idx_all.isin(sample_of_interest)]
            regular_idx = idx_all[~idx_all.isin(sample_of_interest)]
            cats = metadata.loc[regular_idx, meta_col].dropna()
            counts = cats.value_counts()
            count_map = {cat: int(counts.get(cat, 0)) for cat in all_cats}
            freq_map = {cat: (count_map[cat] / total_per_cat[cat] if cat in total_per_cat and total_per_cat[cat] > 0 else 0.0) for cat in all_cats}

            if len(cats) > 0:
                if mode == "frequency":
                    dominant_cat = max(freq_map, key=freq_map.get)
                else:
                    dominant_cat = max(count_map, key=count_map.get)
                type_val = ("exclusive" if sum(1 for v in count_map.values() if v > 0) == 1 else "shared")

            elif len(soi_idx) > 0:
                dominant_cat = ("sample_of_interest")
                if len(soi_idx) > 1:
                    type_val = ("exclusive " "(group sample_of_interest)")
                else:
                    type_val = ("exclusive " "(this sample_of_interest)")

            else:
                dominant_cat = "-"
                type_val = "unclassified"
            counts_str = "; ".join([f"{cat}:{count_map[cat]}" for cat in sorted(all_cats)])
            freqs_str = "; ".join([f"{cat}:{freq_map[cat]:.4f}" for cat in sorted(all_cats)])
            rows.append([locus, allele, dominant_cat, type_val, counts_str, freqs_str])

    df = pd.DataFrame(rows, columns=["Locus", "Allele", "Dominant_Category", "Type", "Counts_per_category", "Frequencies_per_category"])
    out = (f"{output_prefix}_allele_dominance.tsv")
    df.to_csv(out, sep="\t", index=False)
    print_log(f"\t[OK] Allele dominance table saved to {out}", log)

def main():
    # define arguments
    parser = argparse.ArgumentParser(description=("Fast allelic heatmap ordered by tree"))
    parser.add_argument("-a", "--allelic_matrix", required=True, help="Allelic matrix in TSV format.")
    parser.add_argument("-m", "--metadata", required=True, help="Metadata table in TSV format.")
    parser.add_argument("-c", "--colors", required=False, help="(Optional) TSV file indicating the colors to use in the heatmap. This file should contain the following columns in the header: 'column' (indicating the name of \
                        the metadata column to use for categories), 'category' (indicating the different categories observed in the specified metadata column) and color (indicating the color code to use for the category). \
                        In case you specify '--sample_of_interest', you must also indicate the color you want to use for alleles exclusive of samples of interest. If this argument is omitted, colors are auto-assigned and \
                        saved to <prefix>_colors.tsv.")
    parser.add_argument("-hcol", "--header", required=True, help="Metadata column to use for categories.")
    parser.add_argument("-t", "--tree", required=True, help="Newick tree file.")
    parser.add_argument("-o", "--output_prefix", dest="output", required=True, help="Output prefix.")
    parser.add_argument("-mode", "--mode", choices=["frequency", "count"], default="frequency", help="Dominant category mode: 'frequency' (default) or 'count'")
    parser.add_argument("-ncat", "--ncat", type=int, default=None, help="Optional: if an allele is present in at least N categories, color it gray (nearly conserved).")
    parser.add_argument("-ncat_dominant", "--ncat_dominant", type=int, default=None, help="Optional: if an allele is dominant in at least N categories, color it gray (nearly conserved)")
    parser.add_argument("--sample_of_interest", required=False, help="[OPTIONAL] TSV or TXT file with the list of samples of interest. Sample IDs must be in the first column and no headers should be provided.")
    parser.add_argument("--tree-labels", type=float, default=None, metavar="SIZE", help="[OPTIONAL] Show sample names at the tree tips. SIZE is the font size in points (e.g. 6).")
    parser.add_argument("--tree-scale", type=float, default=None, metavar="SIZE", help="[OPTIONAL] Show the horizontal tree branch-length scale " "at the bottom. SIZE is the tick-label font size in points " "(e.g. 8).")

    args = parser.parse_args()

    log_name = args.output + ".log"
    log = open(log_name, "w+")

    print_log("\n******************** running allelic_heatmap.py ********************\n", log)
    print_log("version " + str(version) + " last updated on " + str(last_updated) + "\n", log)
    print_log(" ".join(sys.argv), log)

    start = datetime.datetime.now()
    print_log("start: " + str(start), log)

    allelic_matrix = pd.read_csv(args.allelic_matrix, sep="\t", index_col=0, dtype=str)
    metadata = pd.read_csv(args.metadata, sep="\t", index_col=0, dtype=str)
    sample_order = read_tree_order(args.tree)

    if args.sample_of_interest:
        print_log("\nSample of interest provided...", log)
        sample_of_interest = (read_sample_of_interest(args.sample_of_interest))
        sample_of_interest_present = (sample_of_interest .intersection(set(allelic_matrix.index)))
        missing_soi = (sample_of_interest - set(allelic_matrix.index))
        if missing_soi:
            print_log("\t[WARNING] The following samples " "of interest were not found in " "the allelic matrix and will be ignored:", log)
            for sample in sorted(missing_soi):
                print_log(f"\t\t  {sample}", log)
        sample_of_interest = (sample_of_interest_present)
        print_log(f"\t[INFO] Samples of interest: " f"{len(sample_of_interest)}", log)
    else:
        sample_of_interest = set()

    sample_interest_color = SAMPLE_INTEREST_COLOR

    if args.colors:
        print_log("\nColors provided...", log)

        colors_df = pd.read_csv(args.colors, sep="\t", dtype=str)
        csub = colors_df[colors_df["column"] == args.header]

        if csub.empty:
            raise ValueError(
                f"No colors found in {args.colors} for column '{args.header}'."
            )

        # Read the special sample_of_interest colour separately.
        soi_rows = csub[csub["category"] == "sample_of_interest"]

        if not soi_rows.empty:
            sample_interest_color = soi_rows.iloc[0]["color"]
        elif sample_of_interest:
            raise ValueError(
                "Your colours file must contain a row with "
                "category='sample_of_interest' for the selected metadata column."
            )

        # sample_of_interest is NOT a metadata category.
        # Therefore, exclude it from the ordinary category colour map.
        color_map = {
            row["category"]: row["color"]
            for _, row in csub.iterrows()
            if row["category"] != "sample_of_interest"
        }

    else:
        regular_metadata = metadata.loc[
            ~metadata.index.isin(sample_of_interest)
        ]

        cats = regular_metadata[args.header].dropna().unique()

        color_map = generate_auto_colors(cats)

        sample_interest_color = SAMPLE_INTEREST_COLOR

        save_auto_colors_tsv(
            log,
            color_map,
            args.output,
            args.header,
            sample_interest_color
        )

    (allele_to_color, exclusives, shared, sample_interest_exclusives, allele_to_categories, conserved_loci) = assign_colors(allelic_matrix, metadata, args.header, color_map, args.mode, args.ncat, args.ncat_dominant, sample_of_interest, sample_interest_color)
    print_log("\nGenerating the heatmap...", log)
    build_heatmap(log, allelic_matrix, allele_to_color, allele_to_categories, sample_order, color_map, args.header, args.output, sample_of_interest, sample_interest_exclusives, args.tree, args.tree_labels, args.tree_scale, sample_interest_color)
    print_log("\nSaving summary stats...", log)
    save_summary(log, exclusives, shared, sample_interest_exclusives, allele_to_categories, args.output)
    print_log("\nBuilding detailed matrix...", log)
    build_detailed_matrix(log, allelic_matrix, sample_order, allele_to_categories, conserved_loci, args.output, sample_interest_exclusives)
    print_log("\nBuilding dominance table...", log)
    build_dominance_table(log, allelic_matrix, metadata, args.header, args.mode, args.output, sample_of_interest)

    end = datetime.datetime.now()
	
    print_log("\nend: " + str(end), log)
    print_log("\nDone!", log)
    print_log("\n************************************************************", log)
    
if __name__ == "__main__":
    main()