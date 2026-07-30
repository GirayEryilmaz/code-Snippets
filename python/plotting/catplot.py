import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc


def categorical_dotplot(
    adata,
    value,
    splitby,
    groupby,
    layer=None,
    expression_threshold=0,
    figsize=(8, 5),
    min_dot_size=20,
    max_dot_size=300,
    dodge_width=0.28,
    ax=None,
):
    """
    Plot group-level mean values with:
      - x-axis: groupby
      - y-axis: mean expression/value
      - color: splitby
      - dot size: percentage of cells above expression_threshold
    """

    df = sc.get.obs_df(
        adata,
        keys=[value, splitby, groupby],
        layer=layer,
    ).rename(columns={value: "_value"})

    stats = (
        df.groupby([groupby, splitby], observed=True)["_value"]
        .agg(
            mean="mean",
            pct=lambda x: 100 * (x > expression_threshold).mean(),
        )
        .reset_index()
    )

    if isinstance(df[groupby].dtype, pd.CategoricalDtype):
        group_order = list(df[groupby].cat.categories)
    else:
        group_order = list(pd.unique(df[groupby]))

    if isinstance(df[splitby].dtype, pd.CategoricalDtype):
        split_order = list(df[splitby].cat.categories)
    else:
        split_order = list(pd.unique(df[splitby]))

    group_positions = dict(zip(group_order, range(len(group_order))))

    # Stable, partial horizontal separation of splitby categories
    if len(split_order) == 1:
        offsets = np.array([0.0])
    else:
        offsets = np.linspace(
            -dodge_width / 2,
            dodge_width / 2,
            len(split_order),
        )

    colors = plt.get_cmap("tab10").colors

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    for i, split in enumerate(split_order):
        sub = stats.loc[stats[splitby] == split]

        x = (
            sub[groupby]
            .map(group_positions)
            .astype(float)
            + offsets[i]
        )

        sizes = (
            min_dot_size
            + sub["pct"] / 100 * (max_dot_size - min_dot_size)
        )

        ax.scatter(
            x,
            sub["mean"],
            s=sizes,
            color=colors[i % len(colors)],
            label=str(split),
            alpha=0.85,
            edgecolor="black",
            linewidth=0.4,
        )

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, rotation=45, ha="right")
    ax.set_xlim(-0.5, len(group_order) - 0.5)

    ax.set_xlabel(groupby)
    ax.set_ylabel(f"Mean {value}")

    split_legend = ax.legend(
        title=splitby,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(split_legend)

    pct_values = [25, 50, 75, 100]

    size_handles = [
        ax.scatter(
            [],
            [],
            s=(
                min_dot_size
                + pct / 100 * (max_dot_size - min_dot_size)
            ),
            facecolor="lightgray",
            edgecolor="black",
            linewidth=0.4,
            label=f"{pct}%",
        )
        for pct in pct_values
    ]

    ax.legend(
        handles=size_handles,
        title="% expressing",
        bbox_to_anchor=(1.02, 0),
        loc="lower left",
        frameon=False,
    )

    ax.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(right=0.78, bottom=0.2)

    return fig, ax, stats
