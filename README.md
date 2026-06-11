# Genómica

Proyecto final de Genómica Computacional basado en la serie pública
`GSE157103` de GEO/NCBI. El objetivo es evaluar si perfiles de expresión génica
en sangre permiten identificar subgrupos de pacientes con COVID-19 relacionados
con el estado UCI/no UCI.

## Integrantes

- Rivera Morales David
- Lopez Bernal Yeimi Lizet
- Hidalgo Carrillo Amir Gilberto
- Badager Estrada Aaron Omar

## Resultado principal

El mejor método fue agrupamiento jerárquico con enlace Ward y `k=2`.

| Métrica | Valor |
|---|---:|
| ARI | 0.263 |
| NMI | 0.207 |
| Silueta | 0.176 |
| Exactitud binaria aproximada | 0.76 |

La interpretación correcta es que existe una señal molecular parcial de
severidad. No es un clasificador clínico perfecto, pero el grupo enriquecido en
pacientes UCI también muestra genes y procesos compatibles con inflamación,
respuesta inmune innata y degranulación de neutrófilos.

## Estructura

| Ruta | Contenido |
|---|---|
| `scripts/` | Pipeline de análisis y generación de PDFs. |
| `config/` | Parámetros principales del proyecto. |
| `data/raw/README_datos.md` | Descripción de los datos y descarga reproducible. |
| `results/tables/` | Tablas principales de resultados. |
| `results/figures/` | Figuras usadas en reporte y presentación. |
| `report/` | Reporte final en Markdown, LaTeX y PDF. |
| `presentation/` | Presentación PDF, fuente Markdown y guion. |
| `docs/fuentes_verificadas.md` | Fuentes revisadas para datos, métodos y referencias. |
| `tests/` | Pruebas básicas del pipeline. |

## Reproducción

```bash
make init
make analisis
make pdfs
make check
```

Si los archivos de GEO no existen en `data/raw/`, el script
`scripts/run_analysis.py` intenta descargarlos desde las URLs oficiales.

Para compilar el reporte PDF desde LaTeX se requiere `pdflatex`.

## Entregables

- `report/reporte_final.pdf`
- `presentation/presentacion.pdf`
- `presentation/guion_presentacion.md`

## Fuentes principales

- GEO/NCBI: `GSE157103`
- Overmyer et al., *Cell Systems*, 2021.
- Pedregosa et al., *Journal of Machine Learning Research*, 2011.
- McInnes et al., *Journal of Open Source Software*, 2018.
- Raudvere et al., *Nucleic Acids Research*, 2019.
- Benjamini y Hochberg, *Journal of the Royal Statistical Society: Series B*,
  1995.
