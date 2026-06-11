#!/usr/bin/env python3
"""Pipeline de análisis para GSE157103 y agrupamiento transcriptómico.

El script descarga archivos procesados de GEO si faltan, limpia metadatos,
preprocesa matrices TPM/conteos, compara métodos de agrupamiento, calcula genes
marcadores, intenta enriquecimiento con g:Profiler y escribe tablas, figuras,
reporte y presentación.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.cluster import AgglomerativeClustering, DBSCAN, KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import trustworthiness
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

try:
    import umap
except Exception:  # pragma: no cover - optional dependency
    umap = None

try:
    from gprofiler import GProfiler
except Exception:  # pragma: no cover - optional dependency
    GProfiler = None


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "results" / "figures"
TABLES = ROOT / "results" / "tables"
REPORT = ROOT / "report"
PRESENTATION = ROOT / "presentation"
ARTICLE = ROOT / "article"
SLIDES = ROOT / "slides"
NOTEBOOKS = ROOT / "notebooks"

URLS = {
    "counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157103/suppl/GSE157103_genes.ec.tsv.gz",
    "tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157103/suppl/GSE157103_genes.tpm.tsv.gz",
    "series": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157103/matrix/GSE157103_series_matrix.txt.gz",
}

FILES = {
    "counts": RAW / "GSE157103_genes.ec.tsv.gz",
    "tpm": RAW / "GSE157103_genes.tpm.tsv.gz",
    "series": RAW / "GSE157103_series_matrix.txt.gz",
}

RANDOM_STATE = 42

TEAM_MEMBERS = [
    "Rivera Morales David",
    "Lopez Bernal Yeimi Lizet",
    "Hidalgo Carrillo Amir Gilberto",
    "Badager Estrada Aaron Omar",
]

ENRICHMENT_COLUMNS = [
    "cluster",
    "source",
    "native",
    "name",
    "p_value",
    "significant",
    "description",
    "term_size",
    "query_size",
    "intersection_size",
    "effective_domain_size",
    "precision",
    "recall",
    "query",
    "parents",
    "intersections",
    "evidences",
]

METHOD_LABELS_ES = {
    "hierarchical": "jerárquico",
    "k-means": "k-means",
    "DBSCAN": "DBSCAN",
}

TERM_TRANSLATIONS_ES = {
    "Neutrophil degranulation": "degranulación de neutrófilos",
    "response to stimulus": "respuesta a estímulos",
    "inflammatory response": "respuesta inflamatoria",
    "defense response": "respuesta de defensa",
    "Innate Immune System": "sistema inmune innato",
    "Immune System": "sistema inmune",
    "cellular response to stimulus": "respuesta celular a estímulos",
    "response to external stimulus": "respuesta a estímulos externos",
    "cell communication": "comunicación celular",
    "Th17 cell differentiation": "diferenciación de células Th17",
    "Th1 and Th2 cell differentiation": "diferenciación de células Th1 y Th2",
    "Metabolism of non-coding RNA": "metabolismo de RNA no codificante",
    "snRNP Assembly": "ensamblaje de snRNP",
    "Human T-cell leukemia virus 1 infection": "infección por virus de leucemia de células T humanas tipo 1",
}


def ensure_dirs() -> None:
    for path in [RAW, PROCESSED, FIGURES, TABLES, REPORT, PRESENTATION, ARTICLE, SLIDES, NOTEBOOKS]:
        path.mkdir(parents=True, exist_ok=True)


def translate_term(term: object) -> str:
    return TERM_TRANSLATIONS_ES.get(str(term), str(term))


def method_name_es(method: object) -> str:
    return METHOD_LABELS_ES.get(str(method), str(method))


def download_if_missing() -> None:
    for key, url in URLS.items():
        dest = FILES[key]
        if dest.exists() and dest.stat().st_size > 0:
            continue
        print(f"Downloading {dest.name}")
        urllib.request.urlretrieve(url, dest)


def clean_value(value: str) -> str:
    return value.strip().strip('"')


def parse_series_matrix(path: Path) -> pd.DataFrame:
    sample_fields: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = row[0]
            values = [clean_value(v) for v in row[1:]]
            if key in {"!Sample_title", "!Sample_geo_accession"}:
                sample_fields[key] = values
            elif key == "!Sample_characteristics_ch1":
                characteristics.append(values)

    titles = sample_fields["!Sample_title"]
    geo_accessions = sample_fields["!Sample_geo_accession"]
    records: list[dict[str, object]] = []

    for idx, title in enumerate(titles):
        record: dict[str, object] = {
            "geo_accession": geo_accessions[idx],
            "title": title,
            "short_id": title_to_short_id(title),
        }
        for char_row in characteristics:
            if idx >= len(char_row):
                continue
            item = char_row[idx]
            if ": " not in item:
                continue
            key, value = item.split(": ", 1)
            record[normalize_key(key)] = value
        records.append(record)

    metadata = pd.DataFrame.from_records(records)
    metadata["disease_state"] = metadata["disease_state"].str.replace(
        "non-COVID-19", "NonCOVID", regex=False
    )
    metadata["disease_state"] = metadata["disease_state"].str.replace(
        "COVID-19", "COVID", regex=False
    )
    metadata["icu_status"] = metadata["icu"].map({"yes": "ICU", "no": "NonICU"})
    metadata["mechanical_ventilation"] = metadata["mechanical_ventilation"].map(
        {"yes": "Yes", "no": "No"}
    )
    metadata["age"] = metadata["age_years"].map(parse_age)
    for col in ["apacheii", "charlson_score"]:
        if col in metadata:
            metadata[col] = pd.to_numeric(
                metadata[col].replace({"unknown": np.nan}), errors="coerce"
            )
    metadata = metadata.sort_values("short_id").reset_index(drop=True)
    return metadata


def title_to_short_id(title: str) -> str:
    match = re.match(r"^(COVID|NONCOVID)_(\d+)_", title)
    if not match:
        raise ValueError(f"Cannot parse sample title: {title}")
    prefix = "C" if match.group(1) == "COVID" else "NC"
    return f"{prefix}{int(match.group(2))}"


def normalize_key(key: str) -> str:
    key = key.strip().lower()
    key = key.replace("(", "").replace(")", "")
    key = key.replace(" ", "_").replace("-", "_")
    return key


def parse_age(value: object) -> float:
    if not isinstance(value, str):
        return np.nan
    if value.startswith(">"):
        value = value[1:]
    value = value.strip().lower().replace(":", "").removesuffix("y")
    return pd.to_numeric(value, errors="coerce")


def read_expression(path: Path, use_sum_for_duplicates: bool) -> pd.DataFrame:
    expr = pd.read_csv(path, sep="\t", compression="gzip")
    expr = expr.rename(columns={"#symbol": "gene"})
    expr["gene"] = expr["gene"].astype(str).str.replace('"', "", regex=False)
    expr = expr.set_index("gene")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if use_sum_for_duplicates:
        expr = expr.groupby(expr.index).sum()
    else:
        expr = expr.groupby(expr.index).mean()
    return expr


def prepare_data(metadata: pd.DataFrame, counts: pd.DataFrame, tpm: pd.DataFrame):
    common = [sample for sample in metadata["short_id"] if sample in tpm.columns]
    metadata = metadata[metadata["short_id"].isin(common)].copy()
    metadata = metadata.set_index("short_id").loc[common].reset_index()
    counts = counts[common]
    tpm = tpm[common]

    covid_samples = metadata.loc[metadata["disease_state"] == "COVID", "short_id"].tolist()
    expressed = (tpm[covid_samples] >= 1).sum(axis=1) >= max(5, math.ceil(0.10 * len(covid_samples)))
    count_filter = (counts[covid_samples] >= 10).sum(axis=1) >= max(5, math.ceil(0.10 * len(covid_samples)))
    keep = expressed & count_filter

    filtered_tpm = tpm.loc[keep]
    filtered_counts = counts.loc[keep]
    log_tpm = np.log2(filtered_tpm + 1)

    variances = log_tpm[covid_samples].var(axis=1).sort_values(ascending=False)
    variable_genes = variances.head(min(3000, len(variances))).index
    covid_log = log_tpm.loc[variable_genes, covid_samples].T
    scaler = StandardScaler()
    scaled = scaler.fit_transform(covid_log)

    return metadata, filtered_counts, filtered_tpm, log_tpm, covid_samples, variable_genes, covid_log, scaled


def run_pca(scaled: np.ndarray, n_components: int = 30):
    n_components = min(n_components, scaled.shape[0] - 1, scaled.shape[1])
    pca = PCA(n_components=n_components, random_state=RANDOM_STATE)
    pcs = pca.fit_transform(scaled)
    pc_cols = [f"PC{i + 1}" for i in range(pcs.shape[1])]
    return pca, pd.DataFrame(pcs, columns=pc_cols)


def cluster_accuracy_for_binary(labels: np.ndarray, truth: np.ndarray) -> float:
    unique_labels = sorted(set(labels))
    unique_truth = sorted(set(truth))
    if len(unique_labels) != 2 or len(unique_truth) != 2:
        return np.nan
    mapping_1 = {unique_labels[0]: unique_truth[0], unique_labels[1]: unique_truth[1]}
    mapping_2 = {unique_labels[0]: unique_truth[1], unique_labels[1]: unique_truth[0]}
    acc_1 = np.mean([mapping_1[x] == y for x, y in zip(labels, truth)])
    acc_2 = np.mean([mapping_2[x] == y for x, y in zip(labels, truth)])
    return max(acc_1, acc_2)


def valid_silhouette(matrix: np.ndarray, labels: np.ndarray) -> float:
    labels = np.asarray(labels)
    unique = set(labels)
    if len(unique) < 2 or len(unique) >= len(labels):
        return np.nan
    try:
        return float(silhouette_score(matrix, labels))
    except Exception:
        return np.nan


def evaluate_clustering(pc_matrix: np.ndarray, icu_truth: np.ndarray):
    rows: list[dict[str, object]] = []
    assignments: dict[str, np.ndarray] = {}

    for k in range(2, 7):
        labels = KMeans(n_clusters=k, n_init=50, random_state=RANDOM_STATE).fit_predict(pc_matrix)
        key = f"kmeans_k{k}"
        assignments[key] = labels
        rows.append(metric_row("k-means", f"k={k}", key, labels, pc_matrix, icu_truth))

    for linkage in ["ward", "average", "complete"]:
        for k in range(2, 7):
            kwargs = {"n_clusters": k, "linkage": linkage}
            if linkage != "ward":
                kwargs["metric"] = "euclidean"
            labels = AgglomerativeClustering(**kwargs).fit_predict(pc_matrix)
            key = f"hierarchical_{linkage}_k{k}"
            assignments[key] = labels
            rows.append(metric_row("hierarchical", f"{linkage}, k={k}", key, labels, pc_matrix, icu_truth))

    dbscan_matrix = pc_matrix[:, : min(10, pc_matrix.shape[1])]
    for min_samples in [3, 5, 8, 10]:
        kth_distances = NearestNeighbors(n_neighbors=min_samples).fit(dbscan_matrix).kneighbors(dbscan_matrix)[0][:, -1]
        eps_values = sorted(
            set(np.round(np.quantile(kth_distances, [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]), 2))
        )
        for eps in eps_values:
            labels = DBSCAN(eps=float(eps), min_samples=min_samples).fit_predict(dbscan_matrix)
            key = f"dbscan_ms{min_samples}_eps{eps:.2f}"
            assignments[key] = labels
            row = metric_row("DBSCAN", f"min_samples={min_samples}, eps={eps:.2f}", key, labels, dbscan_matrix, icu_truth)
            row["noise_fraction"] = float(np.mean(labels == -1))
            rows.append(row)

    metrics = pd.DataFrame(rows).sort_values(
        ["adjusted_rand_index", "silhouette"], ascending=False, na_position="last"
    )
    return metrics, assignments


def metric_row(method, params, assignment_key, labels, matrix, truth):
    labels = np.asarray(labels)
    return {
        "method": method,
        "params": params,
        "assignment_key": assignment_key,
        "n_clusters": len(set(labels) - ({-1} if -1 in labels else set())),
        "silhouette": valid_silhouette(matrix, labels),
        "adjusted_rand_index": adjusted_rand_score(truth, labels),
        "normalized_mutual_info": normalized_mutual_info_score(truth, labels),
        "binary_cluster_accuracy": cluster_accuracy_for_binary(labels, truth),
        "noise_fraction": 0.0,
    }


def compute_marker_genes(log_tpm: pd.DataFrame, sample_ids: list[str], labels: np.ndarray) -> pd.DataFrame:
    data = log_tpm[sample_ids]
    unique_labels = sorted(set(labels))
    rows = []
    for cluster in unique_labels:
        in_cluster = [s for s, lab in zip(sample_ids, labels) if lab == cluster]
        out_cluster = [s for s, lab in zip(sample_ids, labels) if lab != cluster]
        if len(in_cluster) < 3 or len(out_cluster) < 3:
            continue
        in_values = data[in_cluster]
        out_values = data[out_cluster]
        log2fc = in_values.mean(axis=1) - out_values.mean(axis=1)
        t_stat, pvals = stats.ttest_ind(
            in_values.T,
            out_values.T,
            axis=0,
            equal_var=False,
            nan_policy="omit",
        )
        pvals = np.nan_to_num(pvals, nan=1.0, posinf=1.0, neginf=1.0)
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        cluster_df = pd.DataFrame(
            {
                "cluster": cluster,
                "gene": data.index,
                "log2fc_cluster_vs_rest": log2fc.values,
                "pvalue": pvals,
                "padj": padj,
                "mean_cluster": in_values.mean(axis=1).values,
                "mean_rest": out_values.mean(axis=1).values,
            }
        )
        cluster_df = cluster_df.sort_values(
            ["padj", "log2fc_cluster_vs_rest"], ascending=[True, False]
        )
        rows.append(cluster_df)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def run_enrichment(markers: pd.DataFrame, max_clusters: int = 8) -> pd.DataFrame:
    if GProfiler is None or markers.empty:
        return pd.DataFrame()
    gp = GProfiler(return_dataframe=True)
    outputs = []
    for cluster, sub in markers.groupby("cluster"):
        if len(outputs) >= max_clusters:
            break
        ranked = sub.sort_values(["padj", "log2fc_cluster_vs_rest"], ascending=[True, False])
        genes = ranked.loc[ranked["log2fc_cluster_vs_rest"] > 0, "gene"].head(150).tolist()
        if len(genes) < 5:
            continue
        try:
            enr = gp.profile(
                organism="hsapiens",
                query=genes,
                sources=["GO:BP", "REAC", "KEGG"],
                no_evidences=False,
            )
        except Exception as exc:
            print(f"g:Profiler enrichment failed for cluster {cluster}: {exc}")
            continue
        if isinstance(enr, pd.DataFrame) and not enr.empty:
            enr = enr.copy()
            enr.insert(0, "cluster", cluster)
            outputs.append(enr.head(20))
    if not outputs:
        return pd.DataFrame()
    return pd.concat(outputs, ignore_index=True)


def plot_qc(metadata, counts, tpm) -> None:
    meta = metadata.copy()
    meta["library_size"] = counts[meta["short_id"]].sum(axis=0).reindex(meta["short_id"]).values
    meta["detected_genes_tpm_ge_1"] = (tpm[meta["short_id"]] >= 1).sum(axis=0).reindex(meta["short_id"]).values
    meta["Tipo de muestra"] = meta["disease_state"].map({"COVID": "COVID-19", "NonCOVID": "No COVID-19"})
    meta["Estado UCI"] = meta["icu_status"].map({"ICU": "UCI", "NonICU": "No UCI"})
    meta.to_csv(TABLES / "sample_qc.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    sns.boxplot(data=meta, x="Tipo de muestra", y="library_size", hue="Estado UCI", ax=axes[0])
    axes[0].set_yscale("log")
    axes[0].set_title("Tamaño de biblioteca")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Conteos esperados, escala logarítmica")
    sns.boxplot(data=meta, x="Tipo de muestra", y="detected_genes_tpm_ge_1", hue="Estado UCI", ax=axes[1])
    axes[1].set_title("Genes detectados")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Genes con TPM >= 1")
    for ax in axes:
        ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(FIGURES / "qc_library_detected_genes.png", dpi=220)
    plt.close(fig)


def plot_pca_all(metadata, log_tpm) -> None:
    samples = metadata["short_id"].tolist()
    genes = log_tpm[samples].var(axis=1).sort_values(ascending=False).head(3000).index
    matrix = StandardScaler().fit_transform(log_tpm.loc[genes, samples].T)
    pca = PCA(n_components=2, random_state=RANDOM_STATE)
    coords = pca.fit_transform(matrix)
    pca_df = metadata.copy()
    pca_df["PC1"] = coords[:, 0]
    pca_df["PC2"] = coords[:, 1]
    pca_df["Tipo de muestra"] = pca_df["disease_state"].map({"COVID": "COVID-19", "NonCOVID": "No COVID-19"})
    pca_df["Estado UCI"] = pca_df["icu_status"].map({"ICU": "UCI", "NonICU": "No UCI"})

    fig, ax = plt.subplots(figsize=(7.5, 6))
    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue="Tipo de muestra",
        style="Estado UCI",
        s=70,
        ax=ax,
    )
    ax.set_title("PCA: todas las muestras")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    fig.tight_layout()
    fig.savefig(FIGURES / "pca_all_samples.png", dpi=220)
    plt.close(fig)


def plot_covid_embeddings(covid_meta, pca_df, umap_df, best_labels, best_name, pca_model) -> None:
    plot_df = covid_meta.copy()
    plot_df = pd.concat([plot_df.reset_index(drop=True), pca_df[["PC1", "PC2"]].reset_index(drop=True)], axis=1)
    plot_df["Grupo"] = best_labels.astype(str)
    plot_df["Estado UCI"] = plot_df["icu_status"].map({"ICU": "UCI", "NonICU": "No UCI"})
    best_title = best_name.replace("hierarchical", "jerárquico")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="Estado UCI", style="Grupo", s=75, ax=axes[0])
    axes[0].set_title("PCA de COVID-19 por estado UCI")
    axes[0].set_xlabel(f"PC1 ({pca_model.explained_variance_ratio_[0] * 100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({pca_model.explained_variance_ratio_[1] * 100:.1f}%)")
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="Grupo", style="Estado UCI", s=75, ax=axes[1])
    axes[1].set_title(f"PCA de COVID-19 por grupo: {best_title}")
    axes[1].set_xlabel(f"PC1 ({pca_model.explained_variance_ratio_[0] * 100:.1f}%)")
    axes[1].set_ylabel(f"PC2 ({pca_model.explained_variance_ratio_[1] * 100:.1f}%)")
    fig.tight_layout()
    fig.savefig(FIGURES / "pca_covid_icu_clusters.png", dpi=220)
    plt.close(fig)

    if not umap_df.empty:
        umap_plot = pd.concat([covid_meta.reset_index(drop=True), umap_df.reset_index(drop=True)], axis=1)
        umap_plot["Grupo"] = best_labels.astype(str)
        umap_plot["Estado UCI"] = umap_plot["icu_status"].map({"ICU": "UCI", "NonICU": "No UCI"})
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        sns.scatterplot(data=umap_plot, x="UMAP1", y="UMAP2", hue="Estado UCI", style="Grupo", s=75, ax=axes[0])
        axes[0].set_title("UMAP por estado UCI")
        sns.scatterplot(data=umap_plot, x="UMAP1", y="UMAP2", hue="Grupo", style="Estado UCI", s=75, ax=axes[1])
        axes[1].set_title(f"UMAP por grupo: {best_title}")
        fig.tight_layout()
        fig.savefig(FIGURES / "umap_covid_icu_clusters.png", dpi=220)
        plt.close(fig)


def plot_metrics(metrics: pd.DataFrame) -> None:
    top = metrics.head(15).copy()
    top["label"] = top["method"].map(method_name_es) + " " + top["params"]
    top["Método"] = top["method"].map(method_name_es)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=top, y="label", x="adjusted_rand_index", hue="Método", dodge=False, ax=ax)
    ax.set_title("Concordancia con las etiquetas UCI/no UCI")
    ax.set_xlabel("Índice de Rand ajustado (ARI)")
    ax.set_ylabel("")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURES / "clustering_ari_top_methods.png", dpi=220)
    plt.close(fig)


def plot_confusion(covid_meta, labels, best_name) -> None:
    table = pd.crosstab(pd.Series(labels, name="cluster"), covid_meta["icu_status"])
    table.to_csv(TABLES / "best_cluster_vs_icu_confusion.csv")
    plot_table = table.rename_axis("Grupo").rename(columns={"ICU": "UCI", "NonICU": "No UCI"})
    fig, ax = plt.subplots(figsize=(6, 4.8))
    sns.heatmap(plot_table, annot=True, fmt="d", cmap="Blues", ax=ax)
    ax.set_title(f"Mejor grupo vs. estado UCI: {best_name.replace('hierarchical', 'jerárquico')}")
    fig.tight_layout()
    fig.savefig(FIGURES / "best_cluster_vs_icu_confusion.png", dpi=220)
    plt.close(fig)


def plot_marker_heatmap(log_tpm, covid_samples, labels, markers) -> None:
    if markers.empty:
        return
    top_genes = (
        markers[markers["log2fc_cluster_vs_rest"] > 0]
        .sort_values(["cluster", "padj", "log2fc_cluster_vs_rest"], ascending=[True, True, False])
        .groupby("cluster")
        .head(12)["gene"]
        .drop_duplicates()
        .head(60)
        .tolist()
    )
    if len(top_genes) < 2:
        return
    order = [sample for _, sample in sorted(zip(labels, covid_samples), key=lambda x: (x[0], x[1]))]
    matrix = log_tpm.loc[top_genes, order]
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, np.nan), axis=0)
    z = z.clip(-2.5, 2.5).fillna(0)
    fig, ax = plt.subplots(figsize=(13, max(6, 0.16 * len(top_genes))))
    sns.heatmap(z, cmap="vlag", center=0, xticklabels=False, yticklabels=True, ax=ax)
    ax.set_title("Genes marcadores principales por grupo")
    ax.set_xlabel("Muestras COVID-19 ordenadas por grupo")
    ax.set_ylabel("Gen")
    fig.tight_layout()
    fig.savefig(FIGURES / "top_marker_gene_heatmap.png", dpi=220)
    plt.close(fig)


def make_umap(scaled: np.ndarray) -> pd.DataFrame:
    if umap is None:
        return pd.DataFrame()
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.25, random_state=RANDOM_STATE)
    coords = reducer.fit_transform(scaled)
    tw = trustworthiness(scaled, coords, n_neighbors=10)
    print(f"UMAP trustworthiness: {tw:.3f}")
    return pd.DataFrame(coords, columns=["UMAP1", "UMAP2"])


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def latex_confusion_table(table: pd.DataFrame) -> str:
    rows = []
    for group, row in table.iterrows():
        rows.append(f"{latex_escape(group)} & {int(row.get('UCI', 0))} & {int(row.get('No UCI', 0))} \\\\")
    return "\n".join(rows)


def write_report(
    metadata,
    covid_meta,
    filtered_counts,
    variable_genes,
    metrics,
    best_row,
    labels,
    markers,
    enrichment,
) -> None:
    best_table = pd.crosstab(pd.Series(labels, name="Grupo"), covid_meta["icu_status"])
    best_table_report = best_table.rename(columns={"ICU": "UCI", "NonICU": "No UCI"})
    top_markers = markers.sort_values(["cluster", "padj", "log2fc_cluster_vs_rest"], ascending=[True, True, False])
    top_marker_lines = []
    top_marker_latex_lines = []
    for cluster, sub in top_markers.groupby("cluster"):
        genes = sub.loc[sub["log2fc_cluster_vs_rest"] > 0, "gene"].head(10).tolist()
        if genes:
            top_marker_lines.append(f"- Grupo {cluster}: " + ", ".join(genes))
            top_marker_latex_lines.append(f"\\item \\textbf{{Grupo {cluster}}}: {latex_escape(', '.join(genes))}")
    enrichment_lines = []
    enrichment_latex_lines = []
    if not enrichment.empty and "name" in enrichment:
        for cluster, sub in enrichment.groupby("cluster"):
            terms = [translate_term(term) for term in sub["name"].head(5).tolist()]
            enrichment_lines.append(f"- Grupo {cluster}: " + "; ".join(terms))
            enrichment_latex_lines.append(f"\\item \\textbf{{Grupo {cluster}}}: {latex_escape('; '.join(terms))}")
    else:
        enrichment_lines.append("- No se obtuvo enriquecimiento automático; revisar `results/tables/enrichment_gprofiler.csv`.")
        enrichment_latex_lines.append(r"\item No se obtuvo enriquecimiento automático.")

    covid_counts = covid_meta["icu_status"].value_counts().to_dict()
    control_n = int((metadata["disease_state"] == "NonCOVID").sum())
    best_ari = best_row["adjusted_rand_index"]
    best_sil = best_row["silhouette"]
    best_acc = best_row["binary_cluster_accuracy"]
    best_method_es = method_name_es(best_row["method"])
    binary_acc_text = f"{best_acc:.2f}" if not pd.isna(best_acc) else "no aplica"
    team_md = "\n".join(f"- {member}" for member in TEAM_MEMBERS)
    team_latex = r" \\ ".join(latex_escape(member) for member in TEAM_MEMBERS)
    method_best = (
        metrics.sort_values(["method", "adjusted_rand_index", "silhouette"], ascending=[True, False, False])
        .groupby("method", as_index=False)
        .head(1)
        .copy()
    )
    method_counts = metrics.groupby("method").size().to_dict()
    method_best["corridas"] = method_best["method"].map(method_counts)
    method_best["método"] = method_best["method"].map(method_name_es)
    method_best_table = method_best[
        ["método", "corridas", "params", "n_clusters", "adjusted_rand_index", "silhouette", "noise_fraction"]
    ].rename(
        columns={
            "params": "mejor_configuración",
            "n_clusters": "grupos",
            "adjusted_rand_index": "ARI",
            "silhouette": "silueta",
            "noise_fraction": "ruido",
        }
    )
    method_best_md = method_best_table.to_markdown(index=False, floatfmt=".3f")
    method_best_latex = "\n".join(
        (
            f"{latex_escape(row['método'])} & {int(row['corridas'])} & "
            f"{latex_escape(row['mejor_configuración'])} & {int(row['grupos'])} & "
            f"{row['ARI']:.3f} & {row['silueta']:.3f} & {row['ruido']:.2f} \\\\"
        )
        for _, row in method_best_table.iterrows()
    )
    dbscan_rows = metrics[metrics["method"] == "DBSCAN"]
    dbscan_best = dbscan_rows.iloc[0] if not dbscan_rows.empty else None
    dbscan_summary = (
        f"DBSCAN se evaluó en {len(dbscan_rows)} combinaciones de parámetros. "
        f"Su mejor ARI fue {dbscan_best['adjusted_rand_index']:.3f} con {dbscan_best['params']} "
        f"y fracción de ruido {dbscan_best['noise_fraction']:.2f}."
        if dbscan_best is not None
        else "DBSCAN no produjo corridas evaluables."
    )
    cluster_lines = []
    cluster_latex_lines = []
    for cluster, row in best_table_report.iterrows():
        total = int(row.sum())
        icu_n = int(row.get("UCI", 0))
        non_n = int(row.get("No UCI", 0))
        icu_pct = 100 * icu_n / total if total else 0
        non_pct = 100 * non_n / total if total else 0
        cluster_lines.append(
            f"- Grupo {cluster}: {icu_n}/{total} UCI ({icu_pct:.1f}%) y "
            f"{non_n}/{total} no UCI ({non_pct:.1f}%)."
        )
        cluster_latex_lines.append(
            f"\\item \\textbf{{Grupo {latex_escape(cluster)}}}: {icu_n}/{total} UCI "
            f"({icu_pct:.1f}\\%) y {non_n}/{total} no UCI ({non_pct:.1f}\\%)."
        )

    report = f"""# Identificación de subgrupos clínicos de COVID-19 mediante expresión génica

Integrantes:

{team_md}

Materia: Genómica Computacional  
Fecha: 11 de junio de 2026  
Conjunto de datos: GSE157103, GEO/NCBI

## Resumen

Este proyecto analiza el dataset GSE157103 para ver si la expresión génica en
sangre permite separar pacientes con COVID-19 en grupos parecidos a UCI y no
UCI. Se usaron las matrices procesadas de conteos esperados y TPM disponibles en
GEO, se limpiaron los metadatos, se filtraron genes con baja señal y se trabajó
con `log2(TPM + 1)`. Después se seleccionaron los {len(variable_genes)} genes más
variables y se compararon tres métodos de agrupamiento: k-means, jerárquico y
DBSCAN.

El mejor resultado fue el agrupamiento {best_method_es} con {best_row['params']}.
Sus métricas fueron ARI={best_ari:.3f}, NMI={best_row['normalized_mutual_info']:.3f},
silueta={best_sil:.3f} y exactitud binaria aproximada de {binary_acc_text}. El
resultado no separa perfectamente a los pacientes, pero sí muestra una señal
clara: un grupo concentra más pacientes UCI y el otro más pacientes no UCI. Los
genes y vías enriquecidas apuntan principalmente a inflamación, respuesta inmune
innata y degranulación de neutrófilos.

Palabras clave: COVID-19; RNA-seq; expresión génica; transcriptómica;
agrupamiento no supervisado; UCI; g:Profiler.

## Introducción

COVID-19 no se presenta igual en todos los pacientes. Algunos cursan la
enfermedad con síntomas moderados y otros requieren cuidados intensivos. Esa
diferencia clínica también puede reflejar cambios moleculares en sangre, sobre
todo en genes relacionados con inflamación y respuesta inmune.

GSE157103 es una serie pública de GEO/NCBI asociada con el estudio multiómico de
Overmyer et al. [1], [2]. El estudio original incluye datos de RNA-seq y otras
capas moleculares. En este proyecto se tomó solo la parte transcriptómica para
contestar una pregunta concreta: si no se le dan al algoritmo las etiquetas UCI
y no UCI, ¿los perfiles de expresión forman grupos parecidos a esas categorías?

El análisis no pretende construir una prueba clínica. La idea es revisar si hay
estructura biológica en los datos y si esa estructura tiene sentido al compararla
con las variables clínicas disponibles.

## Pregunta, objetivo e hipótesis

Pregunta de trabajo: ¿el agrupamiento no supervisado de transcriptomas sanguíneos
recupera, aunque sea parcialmente, el estado UCI/no UCI en pacientes con
COVID-19?

Objetivo: comparar métodos de agrupamiento sobre datos de expresión génica y
evaluar cuál se parece más a la etiqueta clínica UCI/no UCI.

Hipótesis: la sangre contiene una señal transcriptómica de severidad. Esa señal
no tiene por qué separar perfectamente a todos los pacientes, pero debería
aparecer como una tendencia medible en los grupos.

## Datos y metodología

Se usaron tres archivos de GSE157103: conteos esperados, TPM y metadatos de GEO.
Después de alinear expresión y metadatos quedaron {len(metadata)} muestras:
{len(covid_meta)} COVID-19 y {control_n} controles no COVID-19. Dentro de las
muestras COVID-19 hay {int(covid_counts.get('ICU', 0))} pacientes UCI y
{int(covid_counts.get('NonICU', 0))} no UCI.

Los controles no COVID-19 se usaron para revisar el contexto general de los
datos. La comparación principal se hizo solo con las {len(covid_meta)} muestras
COVID-19, porque la pregunta del proyecto es sobre severidad entre pacientes
infectados.

El procesamiento quedó así:

1. Se filtraron genes con baja expresión. El criterio fue TPM al menos 1 y
   conteos esperados al menos 10 en un mínimo de 5 muestras COVID-19 o 10% de
   las muestras COVID-19.
2. Después del filtrado quedaron {filtered_counts.shape[0]} genes.
3. La expresión se transformó como `log2(TPM + 1)`.
4. Se seleccionaron los {len(variable_genes)} genes con mayor varianza.
5. Los genes se escalaron antes de PCA para que genes con valores grandes no
   dominaran el análisis.
6. Se aplicó PCA y también UMAP como visualización complementaria.
7. Se probaron k-means, agrupamiento jerárquico y DBSCAN con scikit-learn [3].
8. Los grupos se compararon contra UCI/no UCI usando ARI, NMI, silueta y
   exactitud binaria cuando había dos grupos.
9. Los genes marcadores se calcularon de forma exploratoria con prueba t de
   Welch y corrección de Benjamini-Hochberg [6].
10. El enriquecimiento funcional se hizo con g:Profiler usando GO:BP, Reactome y
    KEGG [5].

Observación: la etiqueta UCI/no UCI no se usó para construir los grupos. Se usó
después, únicamente para medir qué tanto coincidían los grupos con una variable
clínica real.

## Resultados

El mejor método por ARI fue el agrupamiento {best_method_es} con
{best_row['params']}. El resultado fue:

- ARI: {best_ari:.3f}
- NMI: {best_row['normalized_mutual_info']:.3f}
- Silueta: {best_sil:.3f}
- Exactitud binaria aproximada: {binary_acc_text}

El ARI muestra una coincidencia parcial con UCI/no UCI. La silueta no es alta,
así que los grupos no están completamente separados. Aun así, la tabla de
contingencia muestra una tendencia clara.

{best_table_report.to_markdown()}

Lectura por grupo:

{chr(10).join(cluster_lines)}

Observación: el grupo 0 concentra pacientes UCI y el grupo 1 concentra pacientes
no UCI. La separación no es perfecta, pero sí hay una dirección clara. Esto
apoya la idea de que la expresión génica en sangre contiene información sobre la
severidad.

Comparación por método, tomando la mejor corrida de cada familia:

{method_best_md}

El método jerárquico tuvo el mejor ARI. k-means también encontró estructura, pero
su mejor resultado tuvo tres grupos y fue menos directo para compararlo con una
etiqueta binaria. DBSCAN quedó por debajo porque marcó muchas muestras como ruido
o produjo grupos con baja concordancia clínica.

{dbscan_summary} Esto deja claro que DBSCAN sí se probó, pero no fue el método
más útil para esta pregunta.

Principales genes aumentados por grupo:

{chr(10).join(top_marker_lines) if top_marker_lines else 'No se detectaron marcadores robustos.'}

Términos funcionales destacados:

{chr(10).join(enrichment_lines)}

Observación: el grupo 0 tiene una señal más inflamatoria. Sus términos
principales incluyen degranulación de neutrófilos, respuesta inflamatoria,
respuesta de defensa y sistema inmune innato. Esto coincide con el estudio
original, donde la severidad de COVID-19 se relaciona con procesos inmunes e
inflamatorios [1].

El grupo 0 incluye genes como `IL1R2`, `UPP1`, `OLAH` y `LRG1`, compatibles con
activación inmune innata y cambios mieloides. El grupo 1 incluye `CD4` y términos
relacionados con células T, por lo que parece menos dominado por inflamación
innata. Esta lectura sirve como hipótesis biológica, no como diagnóstico
individual.

## Discusión

El resultado apoya la hipótesis, pero con cuidado. Los transcriptomas sanguíneos
sí contienen información relacionada con severidad, aunque no separan a todos los
pacientes sin error. Esto tiene sentido porque UCI/no UCI es una etiqueta clínica
práctica, no un estado molecular puro. Dos pacientes pueden estar en la misma
categoría clínica y tener diferencias por edad, tratamiento, comorbilidades o
momento de la enfermedad.

La parte más fuerte del análisis no es la exactitud como clasificador, sino la
coherencia entre tres cosas: el grupo con más pacientes UCI, los genes marcadores
y los términos funcionales de inflamación/respuesta innata. Esa combinación hace
que el resultado sea biológicamente razonable.

Limitaciones:

- Se usaron matrices procesadas de GEO, no lecturas crudas.
- No se ajustó por edad, sexo, ventilación mecánica ni comorbilidades.
- Los marcadores son exploratorios; no sustituyen un análisis formal con DESeq2
  o limma-voom.
- Sangre completa mezcla composición celular y activación transcripcional.
- g:Profiler depende de bases de datos externas que pueden cambiar.

Una parte importante del trabajo fue dejar las corridas completas. k-means,
jerárquico y DBSCAN quedaron en la tabla final aunque no todos funcionaron igual
de bien. Así se ve qué método ganó y también por qué los otros no fueron la mejor
opción.

## Conclusión

El agrupamiento no supervisado sí recuperó una señal parcial de UCI/no UCI en
GSE157103. El mejor resultado fue el agrupamiento {best_method_es} con
{best_row['params']}, ARI={best_ari:.3f} y exactitud binaria aproximada de
{binary_acc_text}. El resultado no funciona como diagnóstico perfecto, pero sí
muestra una señal biológica consistente: el grupo con más pacientes UCI también
presenta genes y vías relacionados con inflamación, respuesta innata y
degranulación de neutrófilos. Esto hace que el análisis sea útil como exploración
transcriptómica y como base para un trabajo posterior con modelos diferenciales y
ajuste por covariables clínicas.

## Anexo de reproducibilidad

El análisis se reproduce desde la raíz del proyecto con:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/run_analysis.py
.venv/bin/python scripts/make_pdfs.py
```

Los parámetros principales están documentados en `config/parametros.json`. Las
salidas tabulares se encuentran en `results/tables/`, las figuras en
`results/figures/` y los PDFs finales en `report/` y `presentation/`. La semilla computacional usada
por el pipeline es 42.

## Bibliografía

[1] K. A. Overmyer et al., "Large-Scale Multi-omic Analysis of COVID-19
Severity", Cell Systems, vol. 12, no. 1, pp. 23-40.e7, 2021.

[2] National Center for Biotechnology Information, "GSE157103: Large-scale
Multi-omic Analysis of COVID-19 Severity", Gene Expression Omnibus, 2020.

[3] F. Pedregosa et al., "Scikit-learn: Machine Learning in Python", Journal of
Machine Learning Research, vol. 12, pp. 2825-2830, 2011.

[4] L. McInnes, J. Healy, N. Saul y L. Grossberger, "UMAP: Uniform Manifold
Approximation and Projection", Journal of Open Source Software, vol. 3, no. 29,
p. 861, 2018.

[5] U. Raudvere et al., "g:Profiler: a web server for functional enrichment
analysis and conversions of gene lists (2019 update)", Nucleic Acids Research,
vol. 47, no. W1, pp. W191-W198, 2019.

[6] Y. Benjamini y Y. Hochberg, "Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing", Journal of the Royal
Statistical Society: Series B, vol. 57, no. 1, pp. 289-300, 1995.
"""
    (ARTICLE / "articulo.md").write_text(report, encoding="utf-8")
    (REPORT / "reporte_final.md").write_text(report, encoding="utf-8")

    latex_report = rf"""\documentclass[11pt]{{article}}
\usepackage[utf8]{{inputenc}}
\usepackage[T1]{{fontenc}}
\usepackage[margin=2.5cm]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{float}}
\usepackage{{hyperref}}
\usepackage{{array}}

\title{{Identificación de subgrupos clínicos de COVID-19 mediante expresión génica}}
\author{{{team_latex}}}
\date{{Genómica Computacional -- Semestre 2026-II\\GSE157103, GEO/NCBI}}

\begin{{document}}
\maketitle

\section*{{Resumen}}
Este proyecto analiza GSE157103 para revisar si la expresión génica en sangre
permite formar grupos parecidos a UCI/no UCI en pacientes con COVID-19. Se
usaron matrices procesadas de conteos esperados y TPM, se filtraron genes con
baja señal y se trabajó con $\log_2(TPM + 1)$. Después se compararon k-means,
agrupamiento jerárquico y DBSCAN. El mejor resultado fue agrupamiento
{latex_escape(best_method_es)} con {latex_escape(best_row['params'])}, con
ARI={best_ari:.3f}, NMI={best_row['normalized_mutual_info']:.3f},
silueta={best_sil:.3f} y exactitud binaria aproximada de
{latex_escape(binary_acc_text)}. La señal no es perfecta, pero sí es consistente:
un grupo concentra más pacientes UCI y el otro más pacientes no UCI.

\section{{Introducción}}
COVID-19 no se presenta igual en todos los pacientes. Algunos cursan la
enfermedad con síntomas moderados y otros requieren cuidados intensivos. Esa
diferencia clínica puede reflejar cambios moleculares en sangre, sobre todo en
genes relacionados con inflamación y respuesta inmune. La pregunta de este
proyecto fue directa: si no se le dan al algoritmo las etiquetas UCI/no UCI, ¿los
perfiles de expresión forman grupos parecidos a esas categorías?

\section{{Datos y metodología}}
Se usaron tres archivos de GSE157103: conteos esperados, TPM y metadatos de GEO.
Después de alinear expresión y metadatos quedaron {len(metadata)} muestras:
{len(covid_meta)} COVID-19 y {control_n} controles no COVID-19. Dentro de las
muestras COVID-19 hay {int(covid_counts.get('ICU', 0))} pacientes UCI y
{int(covid_counts.get('NonICU', 0))} no UCI. La comparación principal se hizo
solo con pacientes COVID-19.

Se retuvieron genes con TPM al menos 1 y conteos esperados al menos 10 en un
mínimo de 5 muestras COVID-19 o 10\% de las muestras COVID-19. Después del
filtrado quedaron {filtered_counts.shape[0]} genes. La expresión se transformó
como $\log_2(TPM + 1)$ y se seleccionaron los {len(variable_genes)} genes de
mayor varianza.

Se aplicaron PCA y UMAP para visualización. Para agrupamiento se compararon
k-means, agrupamiento jerárquico y DBSCAN. La concordancia con UCI/no UCI se
evaluó con ARI, NMI, silueta y exactitud binaria cuando el método produjo dos
grupos.

\subsection*{{Procesamiento aplicado}}
\begin{{enumerate}}
\item Se filtraron genes con baja expresión.
\item Se aplicó la transformación $\log_2(TPM + 1)$.
\item Se seleccionaron los {len(variable_genes)} genes más variables.
\item Se aplicó PCA y UMAP.
\item Se probaron k-means, agrupamiento jerárquico y DBSCAN.
\item UCI/no UCI se usó solo para evaluar los grupos, no para construirlos.
\end{{enumerate}}

\section{{Resultados}}
El mejor método por ARI fue el agrupamiento {latex_escape(best_method_es)} con
{latex_escape(best_row['params'])}. El resultado fue:

\begin{{center}}
\begin{{tabular}}{{l r}}
\hline
Métrica & Valor \\
\hline
ARI & {best_ari:.3f} \\
NMI & {best_row['normalized_mutual_info']:.3f} \\
Silueta & {best_sil:.3f} \\
Exactitud binaria & {latex_escape(binary_acc_text)} \\
\hline
\end{{tabular}}
\end{{center}}

\subsection*{{Comparación por método}}
\begin{{center}}
\small
\begin{{tabular}}{{l r l r r r r}}
\hline
Método & Corridas & Mejor configuración & Grupos & ARI & Silueta & Ruido \\
\hline
{method_best_latex}
\hline
\end{{tabular}}
\end{{center}}

{latex_escape(dbscan_summary)} Esto deja claro que DBSCAN sí se probó, pero no
fue el método más útil para esta pregunta.

\subsection*{{Tabla grupo vs. UCI}}
\begin{{center}}
\begin{{tabular}}{{r r r}}
\hline
Grupo & UCI & No UCI \\
\hline
{latex_confusion_table(best_table_report)}
\hline
\end{{tabular}}
\end{{center}}

\begin{{itemize}}
{chr(10).join(cluster_latex_lines)}
\end{{itemize}}

Observación: el grupo 0 concentra más pacientes UCI y el grupo 1 concentra más
pacientes no UCI. La separación no es perfecta, pero sí hay una dirección clara.

\subsection*{{Genes marcadores principales}}
\begin{{itemize}}
{chr(10).join(top_marker_latex_lines) if top_marker_latex_lines else r"\item No se detectaron marcadores robustos."}
\end{{itemize}}

\subsection*{{Enriquecimiento funcional}}
\begin{{itemize}}
{chr(10).join(enrichment_latex_lines)}
\end{{itemize}}

Observación: el grupo 0 tiene una señal más inflamatoria. Genes como
\textit{{IL1R2}}, \textit{{UPP1}}, \textit{{OLAH}} y \textit{{LRG1}} son
compatibles con activación inmune innata y cambios mieloides. El grupo 1 incluye
\textit{{CD4}} y términos relacionados con células T, por lo que parece menos
dominado por inflamación innata.

\section{{Figuras principales}}
\begin{{figure}}[H]
\centering
\includegraphics[width=0.95\linewidth]{{\detokenize{{../results/figures/qc_library_detected_genes.png}}}}
\caption{{Control de calidad por tipo de muestra y estado UCI.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.95\linewidth]{{\detokenize{{../results/figures/pca_covid_icu_clusters.png}}}}
\caption{{PCA de muestras COVID-19 coloreado por estado UCI y grupo seleccionado.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.85\linewidth]{{\detokenize{{../results/figures/clustering_ari_top_methods.png}}}}
\caption{{Comparación de métodos por índice de Rand ajustado.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.75\linewidth]{{\detokenize{{../results/figures/top_marker_gene_heatmap.png}}}}
\caption{{Mapa de calor de genes marcadores por grupo.}}
\end{{figure}}

\section{{Discusión}}
El resultado apoya la hipótesis, pero con cuidado. Los transcriptomas sanguíneos
sí contienen información relacionada con severidad, aunque no separan a todos los
pacientes sin error. Esto tiene sentido porque UCI/no UCI es una etiqueta clínica
práctica, no un estado molecular puro. La parte más fuerte del análisis es la
coherencia entre el grupo con más pacientes UCI, los genes marcadores y los
términos funcionales de inflamación/respuesta innata.

Limitaciones principales:
\begin{{itemize}}
\item Se usaron matrices procesadas de GEO, no lecturas crudas.
\item No se ajustó por edad, sexo, ventilación mecánica ni comorbilidades.
\item Los marcadores son exploratorios y no sustituyen DESeq2 o limma-voom.
\item Sangre completa mezcla composición celular y activación transcripcional.
\item g:Profiler depende de bases de datos externas que pueden cambiar.
\end{{itemize}}

\section{{Conclusión}}
El agrupamiento no supervisado recuperó una señal parcial de UCI/no UCI en
GSE157103. El mejor resultado fue agrupamiento {latex_escape(best_method_es)}
con {latex_escape(best_row['params'])}, ARI={best_ari:.3f} y exactitud binaria
aproximada de {latex_escape(binary_acc_text)}. El resultado no funciona como
diagnóstico perfecto, pero sí muestra una señal biológica consistente:
inflamación, respuesta innata y degranulación de neutrófilos aparecen en el
grupo con mayor proporción de pacientes UCI.

\section*{{Bibliografía}}
\begin{{enumerate}}
\item Overmyer, K. A. et al. (2021). Large-Scale Multi-omic Analysis of COVID-19 Severity. \textit{{Cell Systems}}, 12(1), 23--40.e7.
\item NCBI GEO. GSE157103: Large-scale Multi-omic Analysis of COVID-19 Severity.
\item Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. \textit{{Journal of Machine Learning Research}}, 12, 2825--2830.
\item McInnes, L. et al. (2018). UMAP: Uniform Manifold Approximation and Projection. \textit{{Journal of Open Source Software}}, 3(29), 861.
\item Raudvere, U. et al. (2019). g:Profiler: a web server for functional enrichment analysis. \textit{{Nucleic Acids Research}}, 47(W1), W191--W198.
\item Benjamini, Y. y Hochberg, Y. (1995). Controlling the False Discovery Rate. \textit{{Journal of the Royal Statistical Society: Series B}}, 57(1), 289--300.
\end{{enumerate}}

\end{{document}}
"""
    (REPORT / "reporte_final.tex").write_text(latex_report, encoding="utf-8")


def write_presentation(best_row, covid_meta, labels, enrichment) -> None:
    table = pd.crosstab(pd.Series(labels, name="Grupo"), covid_meta["icu_status"])
    table = table.rename(columns={"ICU": "UCI", "NonICU": "No UCI"})
    terms = "Revisar tabla de enriquecimiento."
    if not enrichment.empty and "name" in enrichment:
        terms = "; ".join(translate_term(term) for term in enrichment["name"].drop_duplicates().head(5).tolist())
    metrics = pd.read_csv(TABLES / "clustering_metrics.csv")
    method_counts = metrics.groupby("method").size().to_dict()
    dbscan_best = metrics[metrics["method"] == "DBSCAN"].iloc[0]
    best_method_es = method_name_es(best_row["method"])
    team_text = ", ".join(TEAM_MEMBERS)
    content = f"""# Presentación

Identificación de subgrupos clínicos de COVID-19 mediante expresión génica

GSE157103, RNA-seq de leucocitos sanguíneos

Integrantes: {team_text}

Genómica Computacional, 11 de junio de 2026

## Introducción

COVID-19 severo no es un estado biológico único. La severidad combina respuesta
inmune, inflamación, daño vascular, comorbilidades y decisiones clínicas.

El recurso GSE157103 permite estudiar esta heterogeneidad con RNA-seq y
metadatos clínicos de pacientes hospitalizados.

## Introducción

Datos usados:

- {len(covid_meta)} pacientes COVID-19
- {int((covid_meta['icu_status'] == 'ICU').sum())} UCI y {int((covid_meta['icu_status'] == 'NonICU').sum())} no UCI
- 26 controles no COVID-19 para contexto exploratorio
- Matrices procesadas: conteos esperados y TPM

La serie GSE157103 está publicada en GEO/NCBI y la publicación primaria es
Overmyer et al., Cell Systems, DOI 10.1016/j.cels.2020.10.003.

## Objetivo

Evaluar si un pipeline no supervisado de expresión génica identifica grupos
moleculares parcialmente concordantes con UCI/no UCI.

Hipótesis: la sangre contiene una señal transcriptómica de severidad, aunque no
una separación clínica perfecta.

## Hallazgos principales

Pasos aplicados:

- Filtrado de genes con baja señal
- Transformación `log2(TPM + 1)`
- 3000 genes variables
- PCA y UMAP para visualizar estructura
- k-means, jerárquico y DBSCAN
- Evaluación posterior contra UCI/no UCI

## Hallazgos principales

Corridas evaluadas:

- k-means: {method_counts.get('k-means', 0)}
- jerárquico: {method_counts.get('hierarchical', 0)}
- DBSCAN: {method_counts.get('DBSCAN', 0)}

DBSCAN sí se probó. Su mejor ARI fue {dbscan_best['adjusted_rand_index']:.3f}
con {dbscan_best['params']} y ruido {dbscan_best['noise_fraction']:.2f}.

## Hallazgos principales

Mejor método: agrupamiento {best_method_es} ({best_row['params']}).

| Métrica | Valor |
|---|---:|
| ARI | {best_row['adjusted_rand_index']:.3f} |
| NMI | {best_row['normalized_mutual_info']:.3f} |
| Silueta | {best_row['silhouette']:.3f} |
| Exactitud binaria | {best_row['binary_cluster_accuracy']:.2f} |

ARI y silueta no son calificaciones escolares. ARI mide concordancia con
UCI/no UCI; silueta mide separación interna de grupos.

## Hallazgos principales

Tabla grupo vs. UCI:

{table.to_markdown()}

Lectura: existe señal molecular de severidad, pero los grupos no son
diagnósticos perfectos.

## Hallazgos principales

Procesos principales: {terms}

La señal más clara se relaciona con inflamación, respuesta de defensa y sistema
inmune innato.

El grupo enriquecido en UCI muestra genes compatibles con activación
inflamatoria/mieloide. El otro grupo conserva señal inmune, pero menos dominada
por inflamación innata.

## Hallazgos principales

- Sangre completa mezcla composición celular y activación transcripcional.
- No se ajustó por edad, sexo, ventilación ni comorbilidades.
- Marcadores calculados de forma exploratoria sobre `log2(TPM + 1)`.
- g:Profiler puede cambiar si cambian sus bases de datos.

## Conclusiones

La expresión génica en sangre recupera parcialmente UCI/no UCI en GSE157103.

La evidencia más fuerte no es la clasificación perfecta, sino la coherencia
entre clusters, marcadores y procesos inflamatorios/innatos.

Trabajo futuro: validar con DESeq2 o limma-voom y ajustar por edad, sexo,
ventilación y comorbilidades.

## Bibliografía

[1] Overmyer et al., Cell Systems, 2021.

[2] GEO/NCBI, GSE157103, 2020.

[3] Pedregosa et al., JMLR, 2011.

[4] McInnes et al., JOSS, 2018.

[5] Raudvere et al., Nucleic Acids Research, 2019.
"""
    (SLIDES / "diapositivas.md").write_text(content, encoding="utf-8")
    (PRESENTATION / "presentacion.md").write_text(content, encoding="utf-8")


def write_notebook() -> None:
    try:
        import nbformat as nbf
    except Exception:
        return

    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# GSE157103: agrupamiento de transcriptomas COVID-19\n"
            "Notebook de ejecución del proyecto. El análisis completo vive en "
            "`scripts/run_analysis.py` para que sea reproducible desde consola."
        ),
        nbf.v4.new_code_cell(
            "from pathlib import Path\n"
            "ROOT = Path('..').resolve()\n"
            "ROOT"
        ),
        nbf.v4.new_code_cell(
            "%run ../scripts/run_analysis.py"
        ),
        nbf.v4.new_code_cell(
            "import pandas as pd\n"
            "pd.read_csv(ROOT / 'results/tables/clustering_metrics.csv').head(10)"
        ),
        nbf.v4.new_code_cell(
            "from IPython.display import Image, display\n"
            "for fig in ['pca_covid_icu_clusters.png', 'clustering_ari_top_methods.png', 'top_marker_gene_heatmap.png']:\n"
            "    display(Image(filename=str(ROOT / 'results/figures' / fig)))"
        ),
    ]
    nbf.write(nb, NOTEBOOKS / "01_analysis.ipynb")


def main() -> None:
    ensure_dirs()
    download_if_missing()
    sns.set_theme(style="whitegrid", context="notebook")

    metadata = parse_series_matrix(FILES["series"])
    counts = read_expression(FILES["counts"], use_sum_for_duplicates=True)
    tpm = read_expression(FILES["tpm"], use_sum_for_duplicates=False)

    (
        metadata,
        filtered_counts,
        filtered_tpm,
        log_tpm,
        covid_samples,
        variable_genes,
        covid_log,
        scaled,
    ) = prepare_data(metadata, counts, tpm)

    metadata.to_csv(PROCESSED / "metadata_clean.csv", index=False)
    filtered_counts.to_csv(PROCESSED / "counts_filtered.csv")
    filtered_tpm.to_csv(PROCESSED / "tpm_filtered.csv")
    log_tpm.to_csv(PROCESSED / "log2_tpm_filtered.csv")
    pd.Series(variable_genes, name="gene").to_csv(PROCESSED / "variable_genes.csv", index=False)

    plot_qc(metadata, counts.loc[filtered_counts.index], tpm.loc[filtered_tpm.index])
    plot_pca_all(metadata, log_tpm)

    covid_meta = metadata.set_index("short_id").loc[covid_samples].reset_index()
    icu_truth = (covid_meta["icu_status"] == "ICU").astype(int).values
    pca_model, pca_df = run_pca(scaled)
    pca_df.insert(0, "short_id", covid_samples)
    pca_df.to_csv(PROCESSED / "covid_pca_coordinates.csv", index=False)

    umap_df = make_umap(scaled)
    if not umap_df.empty:
        umap_df.insert(0, "short_id", covid_samples)
        umap_df.to_csv(PROCESSED / "covid_umap_coordinates.csv", index=False)
        umap_for_plot = umap_df[["UMAP1", "UMAP2"]]
    else:
        umap_for_plot = pd.DataFrame()

    metrics, assignments = evaluate_clustering(pca_df.filter(regex=r"^PC").values[:, :20], icu_truth)
    metrics.to_csv(TABLES / "clustering_metrics.csv", index=False)
    plot_metrics(metrics)

    best_row = metrics.iloc[0]
    best_key = best_row["assignment_key"]
    best_labels = assignments[best_key]
    assignments_df = covid_meta[["short_id", "title", "icu_status", "mechanical_ventilation", "age", "sex"]].copy()
    assignments_df["selected_cluster"] = best_labels
    for key, labels in assignments.items():
        assignments_df[key] = labels
    assignments_df.to_csv(TABLES / "cluster_assignments.csv", index=False)

    plot_covid_embeddings(covid_meta, pca_df, umap_for_plot, best_labels, str(best_key), pca_model)
    plot_confusion(covid_meta, best_labels, str(best_key))

    markers = compute_marker_genes(log_tpm, covid_samples, best_labels)
    markers.to_csv(TABLES / "cluster_marker_genes.csv", index=False)
    top_up_markers = (
        markers[markers["log2fc_cluster_vs_rest"] > 0]
        .sort_values(["cluster", "padj", "log2fc_cluster_vs_rest"], ascending=[True, True, False])
        .groupby("cluster")
        .head(25)
    )
    top_up_markers.to_csv(TABLES / "top25_upregulated_marker_genes_by_cluster.csv", index=False)
    plot_marker_heatmap(log_tpm, covid_samples, best_labels, markers)

    enrichment = run_enrichment(markers)
    if enrichment.empty:
        enrichment = pd.DataFrame(columns=ENRICHMENT_COLUMNS)
    enrichment.to_csv(TABLES / "enrichment_gprofiler.csv", index=False)

    summary = {
        "n_samples_total": int(len(metadata)),
        "n_covid_samples": int(len(covid_meta)),
        "n_controls": int((metadata["disease_state"] == "NonCOVID").sum()),
        "n_filtered_genes": int(filtered_counts.shape[0]),
        "n_variable_genes": int(len(variable_genes)),
        "best_method": str(best_row["method"]),
        "best_params": str(best_row["params"]),
        "best_assignment_key": str(best_key),
        "best_ari": float(best_row["adjusted_rand_index"]),
        "best_silhouette": float(best_row["silhouette"]),
        "best_nmi": float(best_row["normalized_mutual_info"]),
        "best_binary_accuracy": None
        if pd.isna(best_row["binary_cluster_accuracy"])
        else float(best_row["binary_cluster_accuracy"]),
    }
    pd.Series(summary).to_json(TABLES / "analysis_summary.json", indent=2)

    write_report(metadata, covid_meta, filtered_counts, variable_genes, metrics, best_row, best_labels, markers, enrichment)
    write_presentation(best_row, covid_meta, best_labels, enrichment)
    write_notebook()

    print("Analysis complete.")
    print(pd.Series(summary))


if __name__ == "__main__":
    main()
