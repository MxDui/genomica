# Fuentes verificadas

Última actualización: 2026-06-11

Se revisaron fuentes primarias u oficiales para sostener el reporte y la
presentación.

## Dataset

| Elemento | Fuente | Uso en el proyecto |
|---|---|---|
| GSE157103 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103 | Serie GEO/NCBI usada como fuente de datos. |
| Matrices locales | `data/raw/GSE157103_genes.ec.tsv.gz`, `data/raw/GSE157103_genes.tpm.tsv.gz` | Conteos esperados y TPM usados por el pipeline. |
| Metadatos limpios | `data/processed/metadata_clean.csv` | Confirma 126 muestras alineadas: 100 COVID-19 y 26 controles; en COVID-19 hay 50 UCI y 50 no UCI. |

## Artículo primario

| Referencia | Fuente | Uso en el proyecto |
|---|---|---|
| Overmyer et al. 2021 | https://doi.org/10.1016/j.cels.2020.10.003 | Contexto biológico del estudio multiómico de severidad de COVID-19. |
| PubMed | https://pubmed.ncbi.nlm.nih.gov/33096026/ | Registro biomédico del artículo. |
| PMC | https://pmc.ncbi.nlm.nih.gov/articles/PMC7543711/ | Texto completo para revisar diseño y contexto. |

## Métodos

| Método | Fuente | Uso en el proyecto |
|---|---|---|
| scikit-learn | https://jmlr.org/papers/v12/pedregosa11a.html | PCA, k-means, agrupamiento jerárquico, DBSCAN y métricas. |
| DBSCAN en scikit-learn | https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html | Justifica reportar ruido y sensibilidad a parámetros. |
| UMAP | https://doi.org/10.21105/joss.00861 | Visualización no lineal. |
| g:Profiler | https://doi.org/10.1093/nar/gkz369 | Enriquecimiento funcional. |
| Benjamini-Hochberg | https://doi.org/10.1111/j.2517-6161.1995.tb02031.x | Corrección por pruebas múltiples en genes marcadores. |

## Nota

El estudio original incluye un contexto multiómico amplio. Este proyecto usa el
subconjunto transcriptómico disponible localmente y reporta los números reales
analizados por el pipeline.
