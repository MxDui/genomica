# Datos originales

Este directorio contiene los archivos originales descargados desde GEO/NCBI para
la serie publica `GSE157103`.

## Archivos

| Archivo | Fuente | Uso |
|---|---|---|
| `GSE157103_series_matrix.txt.gz` | GEO series matrix | Metadatos clinicos y de muestra. |
| `GSE157103_genes.ec.tsv.gz` | GEO supplementary file | Conteos esperados por gen. |
| `GSE157103_genes.tpm.tsv.gz` | GEO supplementary file | Expresion normalizada en TPM. |

## Descarga reproducible

Si estos archivos faltan, `scripts/run_analysis.py` intenta descargarlos desde
los enlaces oficiales definidos en la constante `URLS`.

## Condiciones de uso

Los datos pertenecen a sus autores y a GEO/NCBI como repositorio de distribucion.
Este proyecto solo los reutiliza con fines academicos y conserva los nombres de
archivo originales para trazabilidad.

