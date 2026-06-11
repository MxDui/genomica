# Identificación de subgrupos clínicos de COVID-19 mediante expresión génica

Integrantes:

- Rivera Morales David
- Lopez Bernal Yeimi Lizet
- Hidalgo Carrillo Amir Gilberto
- Badager Estrada Aaron Omar

Materia: Genómica Computacional  
Fecha: 11 de junio de 2026  
Conjunto de datos: GSE157103, GEO/NCBI

## Resumen

Este proyecto analiza el dataset GSE157103 para ver si la expresión génica en
sangre permite separar pacientes con COVID-19 en grupos parecidos a UCI y no
UCI. Se usaron las matrices procesadas de conteos esperados y TPM disponibles en
GEO, se limpiaron los metadatos, se filtraron genes con baja señal y se trabajó
con `log2(TPM + 1)`. Después se seleccionaron los 3000 genes más
variables y se compararon tres métodos de agrupamiento: k-means, jerárquico y
DBSCAN.

El mejor resultado fue el agrupamiento jerárquico con ward, k=2.
Sus métricas fueron ARI=0.263, NMI=0.207,
silueta=0.176 y exactitud binaria aproximada de 0.76. El
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
Después de alinear expresión y metadatos quedaron 126 muestras:
100 COVID-19 y 26 controles no COVID-19. Dentro de las
muestras COVID-19 hay 50 pacientes UCI y
50 no UCI.

Los controles no COVID-19 se usaron para revisar el contexto general de los
datos. La comparación principal se hizo solo con las 100 muestras
COVID-19, porque la pregunta del proyecto es sobre severidad entre pacientes
infectados.

El procesamiento quedó así:

1. Se filtraron genes con baja expresión. El criterio fue TPM al menos 1 y
   conteos esperados al menos 10 en un mínimo de 5 muestras COVID-19 o 10% de
   las muestras COVID-19.
2. Después del filtrado quedaron 11924 genes.
3. La expresión se transformó como `log2(TPM + 1)`.
4. Se seleccionaron los 3000 genes con mayor varianza.
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

El mejor método por ARI fue el agrupamiento jerárquico con
ward, k=2. El resultado fue:

- ARI: 0.263
- NMI: 0.207
- Silueta: 0.176
- Exactitud binaria aproximada: 0.76

El ARI muestra una coincidencia parcial con UCI/no UCI. La silueta no es alta,
así que los grupos no están completamente separados. Aun así, la tabla de
contingencia muestra una tendencia clara.

|   Grupo |   UCI |   No UCI |
|--------:|------:|---------:|
|       0 |    40 |       14 |
|       1 |    10 |       36 |

Lectura por grupo:

- Grupo 0: 40/54 UCI (74.1%) y 14/54 no UCI (25.9%).
- Grupo 1: 10/46 UCI (21.7%) y 36/46 no UCI (78.3%).

Observación: el grupo 0 concentra pacientes UCI y el grupo 1 concentra pacientes
no UCI. La separación no es perfecta, pero sí hay una dirección clara. Esto
apoya la idea de que la expresión génica en sangre contiene información sobre la
severidad.

Comparación por método, tomando la mejor corrida de cada familia:

| método     |   corridas | mejor_configuración      |   grupos |   ARI |   silueta |   ruido |
|:-----------|-----------:|:-------------------------|---------:|------:|----------:|--------:|
| DBSCAN     |         44 | min_samples=3, eps=22.71 |        6 | 0.060 |    -0.055 |   0.570 |
| jerárquico |         15 | ward, k=2                |        2 | 0.263 |     0.176 |   0.000 |
| k-means    |          5 | k=3                      |        3 | 0.226 |     0.229 |   0.000 |

El método jerárquico tuvo el mejor ARI. k-means también encontró estructura, pero
su mejor resultado tuvo tres grupos y fue menos directo para compararlo con una
etiqueta binaria. DBSCAN quedó por debajo porque marcó muchas muestras como ruido
o produjo grupos con baja concordancia clínica.

DBSCAN se evaluó en 44 combinaciones de parámetros. Su mejor ARI fue 0.060 con min_samples=3, eps=22.71 y fracción de ruido 0.57. Esto deja claro que DBSCAN sí se probó, pero no fue el método
más útil para esta pregunta.

Principales genes aumentados por grupo:

- Grupo 0: IL1R2, C3orf86, UPP1, OLAH, GRB10, ROPN1L, GMFG, LRG1, FLOT1, RSPH14
- Grupo 1: CD4, TMEM109, AARS2, ST3GAL5, TP53, CHST14, DCP1B, PATZ1, EID2, RSAD1

Términos funcionales destacados:

- Grupo 0: degranulación de neutrófilos; respuesta a estímulos; respuesta inflamatoria; respuesta de defensa; sistema inmune innato
- Grupo 1: diferenciación de células Th17; diferenciación de células Th1 y Th2; metabolismo de RNA no codificante; ensamblaje de snRNP; infección por virus de leucemia de células T humanas tipo 1

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
GSE157103. El mejor resultado fue el agrupamiento jerárquico con
ward, k=2, ARI=0.263 y exactitud binaria aproximada de
0.76. El resultado no funciona como diagnóstico perfecto, pero sí
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
