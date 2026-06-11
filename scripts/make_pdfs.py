#!/usr/bin/env python3
"""Create PDF report and PDF slide deck for the GSE157103 project."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image,
    KeepTogether,
    ListFlowable,
    ListItem,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
ARTICLE_DIR = ROOT / "article"
SLIDES_DIR = ROOT / "slides"
REPORT_DIR = ROOT / "report"
PRESENTATION_DIR = ROOT / "presentation"
DIST_DIR = ROOT / "dist"
FIGURES = ROOT / "results" / "figures"
TABLES = ROOT / "results" / "tables"

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


def translate_term(term: object) -> str:
    return TERM_TRANSLATIONS_ES.get(str(term), str(term))


def read_csv_or_empty(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except EmptyDataError:
        return pd.DataFrame(columns=columns)


def compile_latex_report(tex_path: Path) -> Path:
    for _ in range(2):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", tex_path.name],
            cwd=tex_path.parent,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            check=False,
        )
        if result.returncode != 0:
            tail = result.stdout[-4000:] if result.stdout else "Sin salida de pdflatex."
            raise RuntimeError(f"pdflatex falló al compilar {tex_path}:\n{tail}")
    for suffix in [".aux", ".log", ".out"]:
        aux = tex_path.with_suffix(suffix)
        if aux.exists():
            aux.unlink()
    return tex_path.with_suffix(".pdf")


def styles():
    base = getSampleStyleSheet()
    base.add(
        ParagraphStyle(
            name="TitleCenter",
            parent=base["Title"],
            alignment=TA_CENTER,
            fontName="Helvetica-Bold",
            fontSize=22,
            leading=26,
            spaceAfter=18,
        )
    )
    base.add(
        ParagraphStyle(
            name="ReportHeading1",
            parent=base["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=16,
            leading=20,
            spaceBefore=14,
            spaceAfter=8,
        )
    )
    base.add(
        ParagraphStyle(
            name="ReportHeading2",
            parent=base["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=13,
            leading=16,
            spaceBefore=10,
            spaceAfter=6,
        )
    )
    base.add(
        ParagraphStyle(
            name="BodyJust",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=10.5,
            leading=14,
            alignment=TA_LEFT,
            spaceAfter=6,
        )
    )
    base.add(
        ParagraphStyle(
            name="Small",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=8,
            leading=10,
        )
    )
    base.add(
        ParagraphStyle(
            name="SlideTitle",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=28,
            leading=32,
            alignment=TA_LEFT,
            textColor=colors.HexColor("#1f2937"),
            spaceAfter=16,
        )
    )
    base.add(
        ParagraphStyle(
            name="SlideBody",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=17,
            leading=23,
            spaceAfter=8,
        )
    )
    base.add(
        ParagraphStyle(
            name="SlideSmall",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=11,
            leading=14,
        )
    )
    return base


def esc(text: object) -> str:
    text = str(text)
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`(.*?)`", r"<font name='Courier'>\1</font>", text)
    return text


def md_to_flowables(md_text: str, st) -> list:
    flowables = []
    pending_bullets = []

    def flush_bullets():
        if not pending_bullets:
            return
        flowables.append(
            ListFlowable(
                [ListItem(Paragraph(esc(item), st["BodyJust"])) for item in pending_bullets],
                bulletType="bullet",
                leftIndent=18,
                bulletFontName="Helvetica",
            )
        )
        pending_bullets.clear()

    paragraph_lines = []

    def flush_paragraph():
        if paragraph_lines:
            text = " ".join(line.strip() for line in paragraph_lines).strip()
            flowables.append(Paragraph(esc(text), st["BodyJust"]))
            paragraph_lines.clear()

    lines = md_text.splitlines()
    in_table = False
    table_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("|"):
            flush_paragraph()
            flush_bullets()
            in_table = True
            table_lines.append(stripped)
            continue
        if in_table:
            in_table = False
            add_markdown_table(table_lines, flowables, st)
            table_lines = []

        if not stripped:
            flush_paragraph()
            flush_bullets()
            flowables.append(Spacer(1, 5))
        elif stripped.startswith("# "):
            flush_paragraph()
            flush_bullets()
            title = stripped[2:].strip()
            flowables.append(Paragraph(esc(title), st["TitleCenter"]))
        elif stripped.startswith("## "):
            flush_paragraph()
            flush_bullets()
            flowables.append(Paragraph(esc(stripped[3:].strip()), st["ReportHeading1"]))
        elif re.match(r"^\d+\.\s+", stripped):
            flush_paragraph()
            flush_bullets()
            flowables.append(Paragraph(esc(stripped), st["BodyJust"]))
        elif stripped.startswith("- "):
            flush_paragraph()
            pending_bullets.append(stripped[2:].strip())
        else:
            paragraph_lines.append(stripped)

    flush_paragraph()
    flush_bullets()
    if table_lines:
        add_markdown_table(table_lines, flowables, st)
    return flowables


def add_markdown_table(lines: list[str], flowables: list, st) -> None:
    rows = []
    for idx, line in enumerate(lines):
        if idx == 1 and set(line.replace("|", "").strip()) <= {"-", ":"}:
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        rows.append([Paragraph(esc(cell), st["Small"]) for cell in cells])
    if not rows:
        return
    table = Table(rows, hAlign="LEFT", repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e5e7eb")),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#9ca3af")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    flowables.append(table)
    flowables.append(Spacer(1, 8))


def add_numbered_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#6b7280"))
    canvas.drawRightString(doc.pagesize[0] - 0.55 * inch, 0.35 * inch, f"Página {doc.page}")
    canvas.restoreState()


def scaled_image(path: Path, max_width: float, max_height: float) -> Image:
    img = Image(str(path))
    ratio = min(max_width / img.imageWidth, max_height / img.imageHeight)
    img.drawWidth = img.imageWidth * ratio
    img.drawHeight = img.imageHeight * ratio
    return img


def build_report_pdf() -> Path:
    DIST_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out = DIST_DIR / "articulo_final.pdf"
    tex_source = REPORT_DIR / "reporte_final.tex"
    if tex_source.exists():
        report_pdf = compile_latex_report(tex_source)
        shutil.copyfile(report_pdf, out)
        return out

    st = styles()
    doc = SimpleDocTemplate(
        str(out),
        pagesize=A4,
        rightMargin=0.65 * inch,
        leftMargin=0.65 * inch,
        topMargin=0.6 * inch,
        bottomMargin=0.6 * inch,
    )
    source = ARTICLE_DIR / "articulo.md"
    if not source.exists():
        source = REPORT_DIR / "reporte_final.md"
    md_text = source.read_text(encoding="utf-8")
    flowables = md_to_flowables(md_text, st)
    flowables.append(PageBreak())
    flowables.append(Paragraph("Figuras principales", st["ReportHeading1"]))
    for title, fig in [
        ("Control de calidad", "qc_library_detected_genes.png"),
        ("PCA de muestras COVID-19", "pca_covid_icu_clusters.png"),
        ("Concordancia de métodos", "clustering_ari_top_methods.png"),
        ("Grupo seleccionado vs. UCI", "best_cluster_vs_icu_confusion.png"),
        ("Genes marcadores", "top_marker_gene_heatmap.png"),
    ]:
        flowables.append(KeepTogether([
            Paragraph(esc(title), st["ReportHeading2"]),
            scaled_image(FIGURES / fig, 7.0 * inch, 4.6 * inch),
            Spacer(1, 10),
        ]))
    doc.build(flowables, onFirstPage=add_numbered_footer, onLaterPages=add_numbered_footer)
    shutil.copyfile(out, REPORT_DIR / "reporte_final.pdf")
    return out


def slide_deck_pdf() -> Path:
    st = styles()
    DIST_DIR.mkdir(parents=True, exist_ok=True)
    out = DIST_DIR / "presentacion_final.pdf"
    summary = json.loads((TABLES / "analysis_summary.json").read_text(encoding="utf-8"))
    binary_accuracy = summary.get("best_binary_accuracy")
    binary_accuracy_text = f"{binary_accuracy:.2f}" if binary_accuracy is not None else "no aplica"
    confusion = pd.read_csv(TABLES / "best_cluster_vs_icu_confusion.csv")
    markers = pd.read_csv(TABLES / "top25_upregulated_marker_genes_by_cluster.csv")
    enrichment = read_csv_or_empty(TABLES / "enrichment_gprofiler.csv", ENRICHMENT_COLUMNS)
    metrics = pd.read_csv(TABLES / "clustering_metrics.csv")
    method_counts = metrics.groupby("method").size().to_dict()
    dbscan_best = metrics[metrics["method"] == "DBSCAN"].iloc[0]

    doc = SimpleDocTemplate(
        str(out),
        pagesize=landscape(A4),
        rightMargin=0.55 * inch,
        leftMargin=0.55 * inch,
        topMargin=0.45 * inch,
        bottomMargin=0.45 * inch,
    )
    flowables = []

    def slide(title: str, body: list, image: str | None = None):
        flowables.append(Paragraph(esc(title), st["SlideTitle"]))
        body_items = [Paragraph(esc(item), st["SlideBody"]) for item in body]
        if image:
            content = Table(
                [[body_items, scaled_image(FIGURES / image, 5.0 * inch, 3.7 * inch)]],
                colWidths=[4.4 * inch, 5.4 * inch],
            )
            content.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP")]))
            flowables.append(content)
        else:
            flowables.extend(body_items)
        flowables.append(PageBreak())

    slide(
        "Presentación",
        [
            "Identificación de subgrupos clínicos de COVID-19 mediante expresión génica",
            "Conjunto de datos: GSE157103, RNA-seq de leucocitos de sangre completa",
            "Integrantes: Rivera Morales David, Lopez Bernal Yeimi Lizet, Hidalgo Carrillo Amir Gilberto, Badager Estrada Aaron Omar",
            "Genómica Computacional - 11 de junio de 2026",
        ],
    )
    slide(
        "Introducción",
        [
            "COVID-19 severo no es un estado biológico único: combina respuesta inmune, inflamación, daño vascular y heterogeneidad clínica.",
            f"Muestras totales: {summary['n_samples_total']}; COVID-19: {summary['n_covid_samples']}; controles: {summary['n_controls']}.",
            "La etiqueta clínica principal para evaluar los grupos fue UCI/no UCI.",
        ],
        "qc_library_detected_genes.png",
    )
    slide(
        "Introducción",
        [
            f"Genes retenidos después de filtrado: {summary['n_filtered_genes']}.",
            f"Genes variables usados para PCA y agrupamiento: {summary['n_variable_genes']}.",
            "Transformación principal: log2(TPM + 1), escalada por gen.",
            "GEO/NCBI confirma GSE157103 como serie pública de Homo sapiens con expresión por secuenciación de alto rendimiento.",
            "Los archivos locales de conteos y TPM contienen 126 muestras alineadas.",
        ],
        "pca_all_samples.png",
    )
    slide(
        "Hallazgos principales",
        [
            "Filtrar genes reduce ruido de genes casi apagados.",
            "log2(TPM + 1) comprime valores extremos y ayuda a comparar genes.",
            "PCA resume miles de genes antes de agrupar.",
            "UCI/no UCI se usa después, solo para evaluar concordancia.",
            "La semilla computacional del pipeline es 42.",
        ],
    )
    slide(
        "Objetivo",
        [
            "Evaluar si un pipeline no supervisado de expresión génica identifica grupos moleculares parcialmente concordantes con UCI/no UCI.",
            "Hipótesis: la sangre contiene una señal transcriptómica de severidad, aunque no una separación clínica perfecta.",
            "Se probaron k-means, agrupamiento jerárquico y DBSCAN.",
        ],
        "clustering_ari_top_methods.png",
    )
    slide(
        "Hallazgos principales",
        [
            f"k-means: {method_counts.get('k-means', 0)} corridas con distintos k.",
            f"Jerárquico: {method_counts.get('hierarchical', 0)} corridas con Ward, average y complete.",
            f"DBSCAN: {method_counts.get('DBSCAN', 0)} corridas con búsqueda k-NN de eps y min_samples.",
            f"Mejor DBSCAN: ARI={dbscan_best['adjusted_rand_index']:.3f}, ruido={dbscan_best['noise_fraction']:.2f}.",
            "DBSCAN sí está en la tabla final, pero no fue el método con mayor concordancia clínica.",
        ],
        "clustering_ari_top_methods.png",
    )
    slide(
        "Hallazgos principales",
        [
            f"Mejor método: agrupamiento jerárquico ({summary['best_params']}).",
            f"ARI: {summary['best_ari']:.3f}.",
            f"Silueta: {summary['best_silhouette']:.3f}.",
            f"NMI: {summary['best_nmi']:.3f}.",
            "La recuperación de UCI/no UCI fue parcial, no una clasificación perfecta.",
            "ARI mide concordancia con UCI/no UCI; silueta mide separación interna.",
        ],
        "pca_covid_icu_clusters.png",
    )

    flowables.append(Paragraph("Hallazgos principales", st["SlideTitle"]))
    confusion = confusion.rename(columns={"cluster": "Grupo", "ICU": "UCI", "NonICU": "No UCI"})
    rows = [confusion.columns.tolist()] + confusion.values.tolist()
    table = Table(rows, hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dbeafe")),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#64748b")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 15),
                ("PADDING", (0, 0), (-1, -1), 8),
            ]
        )
    )
    flowables.append(table)
    flowables.append(Spacer(1, 14))
    flowables.append(Paragraph("Grupo 0 concentra 40/54 pacientes UCI; grupo 1 concentra 36/46 pacientes no UCI.", st["SlideBody"]))
    flowables.append(Paragraph("La señal existe, pero hay mezcla clínica relevante.", st["SlideBody"]))
    flowables.append(Spacer(1, 10))
    flowables.append(scaled_image(FIGURES / "best_cluster_vs_icu_confusion.png", 5.0 * inch, 3.0 * inch))
    flowables.append(PageBreak())

    marker_text = []
    for cluster, sub in markers.groupby("cluster"):
        genes = ", ".join(sub["gene"].head(8).tolist())
        marker_text.append(f"Grupo {cluster}: {genes}.")
    slide(
        "Hallazgos principales",
        marker_text + [
            "Grupo 0 muestra genes compatibles con respuesta innata/inflamatoria.",
            "Grupo 1 muestra un perfil más relacionado con células T y regulación celular.",
        ],
        "top_marker_gene_heatmap.png",
    )

    terms = []
    if not enrichment.empty:
        for cluster, sub in enrichment.groupby("cluster"):
            names = "; ".join(translate_term(term) for term in sub["name"].drop_duplicates().head(4).tolist())
            terms.append(f"Grupo {cluster}: {names}.")
    slide(
        "Hallazgos principales",
        (terms[:4] or ["No se obtuvo enriquecimiento automático."])
        + [
            "La señal más clara se relaciona con inflamación, respuesta de defensa y sistema inmune innato.",
            "La interpretación es biológica y exploratoria, no un diagnóstico individual.",
        ],
    )
    slide(
        "Hallazgos principales",
        [
            "Sangre completa mezcla composición celular y activación transcripcional.",
            "No se ajustó por edad, sexo, ventilación mecánica ni comorbilidades.",
            "Los marcadores se calcularon de forma exploratoria sobre log2(TPM + 1).",
            "g:Profiler depende de bases de datos externas que pueden actualizarse.",
        ],
    )
    slide(
        "Conclusiones",
        [
            "Los transcriptomas sanguíneos contienen señal asociada a severidad clínica.",
            f"El mejor agrupamiento recupera parcialmente UCI/no UCI: ARI={summary['best_ari']:.3f} y exactitud binaria aproximada de {binary_accuracy_text}.",
            "La interpretación biológica es más fuerte que una clasificación binaria perfecta: inflamación, respuesta innata y fenotipos inmunes.",
            "Trabajo futuro: validar con DESeq2/limma y ajustar por edad, sexo y comorbilidades.",
        ],
    )
    slide(
        "Bibliografía",
        [
            "Overmyer et al., Cell Systems, 2021.",
            "GEO/NCBI, GSE157103, 2020.",
            "Pedregosa et al., Journal of Machine Learning Research, 2011.",
            "McInnes et al., Journal of Open Source Software, 2018.",
            "Raudvere et al., Nucleic Acids Research, 2019.",
            "Benjamini y Hochberg, JRSS-B, 1995.",
        ],
    )

    doc.build(flowables)
    PRESENTATION_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(out, PRESENTATION_DIR / "presentacion.pdf")
    return out


def main() -> None:
    report = build_report_pdf()
    slides = slide_deck_pdf()
    print(f"Wrote {report}")
    print(f"Wrote {slides}")


if __name__ == "__main__":
    main()
