# Presentación

Identificación de subgrupos clínicos de COVID-19 mediante expresión génica

GSE157103, RNA-seq de leucocitos sanguíneos

Integrantes: Rivera Morales David, Lopez Bernal Yeimi Lizet, Hidalgo Carrillo Amir Gilberto, Badager Estrada Aaron Omar

Genómica Computacional, 11 de junio de 2026

## Introducción

COVID-19 severo no es un estado biológico único. La severidad combina respuesta
inmune, inflamación, daño vascular, comorbilidades y decisiones clínicas.

Esa diferencia clínica puede reflejar cambios moleculares en sangre, sobre todo en genes relacionados con inflamación y respuesta inmune.

El recurso GSE157103 permite estudiar esta heterogeneidad con RNA-seq y
metadatos clínicos de pacientes hospitalizados.

## Datos Verificados

Datos usados:

- 100 pacientes COVID-19
- 50 UCI y 50 no UCI
- 26 controles no COVID-19 para contexto exploratorio
- Matrices procesadas: conteos esperados y TPM

La serie GSE157103 está publicada en GEO/NCBI y la publicación primaria es
Overmyer et al., Cell Systems, DOI 10.1016/j.cels.2020.10.003.

## Objetivo

Evaluar si un pipeline no supervisado de expresión génica identifica grupos
moleculares parcialmente concordantes con UCI/no UCI.

Hipótesis: la sangre contiene una señal transcriptómica de severidad, aunque no
una separación clínica perfecta.

## Pipeline

Pasos aplicados:

- Filtrado de genes con baja señal
- Transformación `log2(TPM + 1)`
- 3000 genes variables
- PCA y UMAP para visualizar estructura
- k-means, jerárquico y DBSCAN
- Evaluación posterior contra UCI/no UCI

## Métodos comparados

Corridas evaluadas:

- k-means: 5
- jerárquico: 15
- DBSCAN: 44

DBSCAN sí se probó. Su mejor ARI fue 0.060
con min_samples=3, eps=22.71 y ruido 0.57.

## Mejor resultado

Mejor método: agrupamiento jerárquico (ward, k=2).

| Métrica | Valor |
|---|---:|
| ARI | 0.263 |
| NMI | 0.207 |
| Silueta | 0.176 |
| Exactitud binaria | 0.76 |

ARI y silueta no son calificaciones escolares. ARI mide concordancia con
UCI/no UCI; silueta mide separación interna de grupos.

## Grupos

Tabla grupo vs. UCI:

|   Grupo |   UCI |   No UCI |
|--------:|------:|---------:|
|       0 |    40 |       14 |
|       1 |    10 |       36 |

Lectura: existe señal molecular de severidad, pero los grupos no son
diagnósticos perfectos.

## Hallazgos principales

Procesos principales: degranulación de neutrófilos; respuesta a estímulos; respuesta inflamatoria; respuesta de defensa; sistema inmune innato

La señal más clara se relaciona con inflamación, respuesta de defensa y sistema
inmune innato.

El grupo enriquecido en UCI muestra genes compatibles con activación
inflamatoria/mieloide. El otro grupo conserva señal inmune, pero menos dominada
por inflamación innata.

## Limitaciones

- Sangre completa mezcla composición celular y activación transcripcional.
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
