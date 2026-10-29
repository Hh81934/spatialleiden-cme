#!/usr/bin/env python3
"""
Original CME (Cross-Method Entropy) script in Python.
"""

import argparse
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def calc_entropy(counts):
    """
    Equivalent to R:
    calc_entropy <- function(u) {
        p <- u[u>0]; p <- p/sum(p); -sum(p*log(p))
    }
    """
    p = counts[counts > 0]
    p = p / p.sum()
    return -np.sum(p * np.log(p))


def align_classes(label_df, ref_col, columns=None):
    """
    Equivalent to R align_classes function.
    Aligns ALL columns simultaneously against reference.
    """
    bcs = label_df.copy()

    # column filtering
    if columns is not None:
        all_cols = bcs.columns.tolist()
        keep_cols = list(set([ref_col] + columns))
        bcs = bcs[keep_cols]

    # convert to numeric-sorted categories
    for col in bcs.columns:
        try:
            unique_vals = sorted(bcs[col].unique(), 
                                  key=lambda x: float(str(x)))
        except ValueError:
            unique_vals = sorted(bcs[col].unique(), key=str)
        bcs[col] = pd.Categorical(bcs[col], categories=unique_vals)

    # get levels
    levs = {col: list(bcs[col].cat.categories) for col in bcs.columns}
    n_levs = {col: len(levs[col]) for col in bcs.columns}

    cols_to_change = [c for c in bcs.columns if c != ref_col]
    ref_labels = bcs[ref_col].astype(str).values
    ref_levels = levs[ref_col]

    for other_col in cols_to_change:

        if n_levs[other_col] != n_levs[ref_col]:
            print(f"{other_col} has {n_levs[other_col]} levels; "
                  f"reference ({ref_col}) has {n_levs[ref_col]} (not modifying)")
            continue

        other_labels = bcs[other_col].astype(str).values
        other_levels = levs[other_col]

        # build contingency table
        cost_matrix = np.zeros((len(ref_levels), len(other_levels)))
        for i, r in enumerate(ref_levels):
            for j, o in enumerate(other_levels):
                cost_matrix[i, j] = np.sum(
                    (ref_labels == str(r)) & (other_labels == str(o))
                )

        # solve Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(-cost_matrix)

        # apply mapping
        new_levels = [None] * len(other_levels)
        for i in range(len(row_ind)):
            new_levels[col_ind[i]] = ref_levels[row_ind[i]]

        level_map = {str(other_levels[j]): str(new_levels[j])
                     for j in range(len(other_levels))
                     if new_levels[j] is not None}

        bcs[other_col] = bcs[other_col].astype(str).map(
            lambda x: level_map.get(x, x)
        )

    # add back filtered columns
    if columns is not None:
        remaining = [c for c in all_cols if c not in bcs.columns]
        if remaining:
            bcs = pd.concat([bcs, label_df[remaining]], axis=1)

    return bcs


def cross_method_entropy(label_df):
    """
    Equivalent to R cross_meth_entropy function.
    Uses align_classes (not align_clusters) to match R behaviour.
    """
    columns = label_df.columns.tolist()
    n_cells = len(label_df)
    entropy_per_ref = []

    for ref_col in columns:

        # align ALL columns at once — matches R:
        # df <- align_classes(d=df, ref=c)
        aligned_df = align_classes(label_df, ref_col)

        # get reference levels — matches R:
        # lv <- levels(df[[c]])
        try:
            ref_unique = np.array(
                sorted(np.unique(aligned_df[ref_col].values),
                       key=lambda x: float(str(x)))
            )
        except ValueError:
            ref_unique = np.sort(np.unique(aligned_df[ref_col].values))

        ref_levels = [str(x) for x in ref_unique]

        # count cluster assignments per cell — matches R:
        # ps <- apply(df, 1, function(u) table(factor(u, levels=lv)))
        # ps <- t(ps)
        counts = np.zeros((n_cells, len(ref_levels)))
        for col in aligned_df.columns:
            col_labels = aligned_df[col].astype(str).values
            for k, cluster in enumerate(ref_levels):
                counts[:, k] += (col_labels == cluster).astype(int)

        # calculate entropy per cell — matches R:
        # ent <- apply(ps, 1, calc_entropy)
        entropy = np.array([
            calc_entropy(counts[i]) for i in range(n_cells)
        ])
        entropy_per_ref.append(entropy)

    # average across all references — matches R:
    # mean_ent <- rowMeans(as.data.frame(cm_ent_list))
    mean_entropy = np.mean(entropy_per_ref, axis=0)
    return mean_entropy


def main():
    parser = argparse.ArgumentParser(
        description="Calculate cross-method entropy at each point."
    )
    parser.add_argument(
        "-i", "--input_file",
        type=str,
        action="append",
        required=True,
        help="Input TSV containing aggregated labels."
    )
    parser.add_argument(
        "-o", "--output_file",
        type=str,
        default="cme_output.tsv",
        help="Output file for spot-wise cross-method entropy."
    )
    parser.add_argument(
        "--BC_ranking",
        type=str,
        default=None,
        help="Optional. Path to TSV or comma-separated BC instances."
    )
    parser.add_argument(
        "-n", "--n_clusters",
        type=int,
        default=None,
        help="Filter input to only use runs with this cluster count."
    )
    parser.add_argument(
        "--n_bcs",
        type=int,
        default=8,
        help="Number of base clusterings to use (default: 8)."
    )
    parser.add_argument(
        "--BC_output",
        type=str,
        default=None,
        help="Save which BC instances were used for CME."
    )
    parser.add_argument(
        "--no_split",
        action="store_true",
        default=False,
        help="Do NOT split by cluster count. Use for control experiment."
    )

    args = parser.parse_args()

    # load all input files
    dfs = []
    for f in args.input_file:
        df = pd.read_csv(f, sep="\t", index_col=0)
        dfs.append(df)
        print(f"Loaded {f}: {df.shape[1]} runs, {df.shape[0]} cells")

    label_df = pd.concat(dfs, axis=1)
    print(f"\nTotal runs: {label_df.shape[1]}")

    if args.BC_ranking is not None:
        import os
        if os.path.exists(args.BC_ranking):
            bc_df = pd.read_csv(args.BC_ranking, sep="\t", index_col=0)
            if args.n_clusters is not None and str(args.n_clusters) in bc_df.columns:
                columns = bc_df[str(args.n_clusters)].dropna().tolist()
            else:
                columns = bc_df.iloc[:, 0].dropna().tolist()
        else:
            columns = args.BC_ranking.split(",")
        columns = columns[:args.n_bcs]
        result_ent = cross_method_entropy(label_df[columns])
        result_df = pd.DataFrame({"entropy": result_ent}, index=label_df.index)

    elif args.no_split:
        columns = label_df.columns[:args.n_bcs].tolist()
        print(f"Running without split ({len(columns)} runs)...")
        result_ent = cross_method_entropy(label_df[columns])
        result_df = pd.DataFrame({"entropy": result_ent}, index=label_df.index)

    else:
        # default — split by cluster count (matches R)
        n_clusters = label_df.apply(lambda x: x.nunique())
        label_lists = {}
        for col, n in n_clusters.items():
            if n not in label_lists:
                label_lists[n] = []
            label_lists[n].append(col)

        print(f"Cluster counts found: {sorted(label_lists.keys())}")

        if args.n_clusters is not None:
            if args.n_clusters in label_lists:
                cols = label_lists[args.n_clusters][:args.n_bcs]
                print(f"Filtering to K={args.n_clusters} ({len(cols)} runs)...")
                result_ent = cross_method_entropy(label_df[cols])
                result_df = pd.DataFrame(
                    {"entropy": result_ent}, index=label_df.index
                )
            else:
                print(f"K={args.n_clusters} not found!")
                return
        else:
            results = {}
            for n, cols in sorted(label_lists.items()):
                cols = cols[:args.n_bcs]
                print(f"Computing CME for K={n} ({len(cols)} runs)...")
                entropy = cross_method_entropy(label_df[cols])
                results[str(n)] = entropy
            result_df = pd.DataFrame(results, index=label_df.index)

    result_df.to_csv(args.output_file, sep="\t")
    print(f"\nSaved: {args.output_file}")
    print(result_df.describe())

    if args.BC_output is not None:
        bc_out = pd.DataFrame({"columns": columns})
        bc_out.to_csv(args.BC_output, sep="\t")
        print(f"Saved BC output: {args.BC_output}")


if __name__ == "__main__":
    main()
