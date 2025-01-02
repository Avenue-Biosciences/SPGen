import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage


def plot_rank_vs_ef(
    df: pd.DataFrame,
    output_path: str,
    ef_col: str,
    ef_rank_col: str,
    protein_name: str = "",
):
    df = df.copy()
    # Different styling for top 20 and 20-100
    df["highlight"] = "Other"
    df.loc[df[ef_rank_col] <= 100, "highlight"] = "Top 100"
    df.loc[df[ef_rank_col] <= 20, "highlight"] = "Top 20"
    df.loc[df["origin"] == "Cytosolic proteins", "highlight"] = "Cytosolic"

    plt.figure(figsize=(10, 6))
    # Plot layers from least important to most important
    for category in ["Other", "Top 100", "Top 20", "Cytosolic"]:
        subset = df[df["highlight"] == category]
        sns.scatterplot(
            x=ef_rank_col,
            y=ef_col,
            data=subset,
            s=12,
            color={
                "Other": "black",
                "Cytosolic": "orange",
                "Top 100": "blue",
                "Top 20": "green",
            }[category],
            edgecolor="none",
            label=category,
        )
    plt.xlabel("Enrichment factor rank")
    plt.ylabel("Enrichment factor")
    handles, labels = plt.gca().get_legend_handles_labels()
    if len(handles) == 4:
        order = [2, 1, 3, 0]
        handles = [handles[idx] for idx in order]
        labels = [labels[idx] for idx in order]

    plt.legend(
        handles,
        labels,
        title="",
        loc="upper right",
    )
    plt.title(f"Enrichment factors for SPs with {protein_name}")
    plt.savefig(output_path, bbox_inches="tight")


def volcano_plot(
    df: pd.DataFrame,
    output_path: str,
    ef_col: str,
    protein_name: str = "",
    top_percent_ef: float = 0.1,
    top_percent_p: float = 0.1,
):
    df = df.copy()
    df["log_p_value"] = -np.log10(df["p_value"])
    df["log_ef"] = np.log2(df[ef_col])

    # Identify top X% in terms of both EF and p-value
    ef_threshold = df["log_ef"].quantile(1 - top_percent_ef)
    p_value_threshold = df["log_p_value"].quantile(1 - top_percent_p)
    df["highlight"] = (df["log_ef"] >= ef_threshold) & (
        df["log_p_value"] >= p_value_threshold
    )

    plt.figure(figsize=(10, 6))

    # Plot non-highlighted points
    sns.scatterplot(
        x="log_ef",
        y="log_p_value",
        data=df[~df["highlight"]],
        s=12,
        color="none",
        edgecolor="black",
        alpha=0.7,
    )

    # Plot highlighted points
    sns.scatterplot(
        x="log_ef",
        y="log_p_value",
        data=df[df["highlight"]],
        s=12,
        color="#D5FF5A",
        edgecolor="black",
        alpha=1,
        label=f"Top {int(top_percent_ef*100)}% EF & {int(top_percent_p*100)}% p-value",
    )

    # Get the maximum absolute value of log_ef to set symmetric limits
    max_abs_log_ef = max(abs(df["log_ef"].min()), abs(df["log_ef"].max()))
    plt.xlim(-max_abs_log_ef * 1.05, max_abs_log_ef * 1.05)

    # Customize y-axis to show original p-values
    min_p_value = df["p_value"].min()
    yticks = [-np.log10(0.05)] + list(range(1, int(-np.log10(min_p_value)) + 1))
    yticklabels = ["0.05"] + [
        f"1e-{i}" for i in range(1, int(-np.log10(min_p_value)) + 1)
    ]
    plt.yticks(yticks, yticklabels)

    # Customize x-axis to show original enrichment factor values
    # Calculate the range of log2 enrichment factors
    min_log_ef = np.floor(df["log_ef"].min())
    max_log_ef = np.ceil(df["log_ef"].max())

    # Create symmetric range
    max_abs_log_ef = int(max(abs(min_log_ef), abs(max_log_ef)))

    # Generate xticks using the largest power of 2
    xticks = [x for x in range(-max_abs_log_ef, max_abs_log_ef + 1)]

    xticks.sort()

    # Create labels
    xticklabels = [f"{int(2**x)}" if x >= 0 else f"1/{int(2**abs(x))}" for x in xticks]

    plt.xticks(xticks, xticklabels)

    plt.axvline(x=0, color="black", linestyle="--")
    plt.xlabel("Enrichment factor")
    plt.ylabel("p-value")
    plt.legend(title="", loc="upper right")
    plt.title(f"Volcano plot for SPs with {protein_name}")
    plt.savefig(output_path, bbox_inches="tight")


def plot_correlation_heatmap(
    df: pd.DataFrame, hf_columns: list[str], lf_columns: list[str], output_path: str
):

    # Correlations between HFs and LFs
    hf_corr = df[hf_columns].corr(method="pearson")
    lf_corr = df[lf_columns].corr(method="pearson")
    corr_matrix = pd.concat([hf_corr, lf_corr])
    plt.figure(figsize=(10, 6))
    sns.heatmap(corr_matrix, annot=True, cmap="coolwarm_r", center=0, vmin=-1, vmax=1)
    plt.title("Correlation heatmap between HFs and LFs")
    plt.savefig(output_path, bbox_inches="tight")


def plot_rank_pairplot(df: pd.DataFrame, rank_cols: list[str], output_path: str):
    # Change output path to PDF
    pdf_path = output_path.replace(".png", ".pdf")

    with PdfPages(pdf_path) as pdf:
        # First pairplot - all data
        g = sns.pairplot(
            df,
            vars=rank_cols,
            corner=True,
            markers="o",
            plot_kws={"alpha": 0.5, "s": 10},
            diag_kind=None,
        )
        g.figure.suptitle(
            "Ranks from different metrics plotted against each other", y=1.02
        )
        pdf.savefig(bbox_inches="tight")
        plt.close()

        # Second pairplot - zoomed to 1-100 range
        g = sns.pairplot(
            df,
            vars=rank_cols,
            corner=True,
            markers="o",
            plot_kws={"alpha": 1, "s": 10},
            diag_kind=None,
        )
        g.figure.suptitle("Zoomed to 1-100 ranks", y=1.02)
        # Set limits for all axes and add dashed lines
        for ax in g.axes.flat:
            if ax is not None:  # Some axes might be None due to corner=True
                ax.set_xlim(0, 100)
                ax.set_ylim(0, 100)
                # Add horizontal and vertical lines that span the entire plot
                ax.axvline(
                    x=20, color="grey", alpha=0.5, linestyle="--"
                )  # vertical line
                ax.axhline(
                    y=20, color="grey", alpha=0.5, linestyle="--"
                )  # horizontal line

        pdf.savefig(bbox_inches="tight")
        plt.close()
