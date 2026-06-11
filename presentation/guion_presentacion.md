# Guion para la presentación

Duración sugerida: 8 a 10 minutos.  
Tono: claro, directo y natural. No conviene leer palabra por palabra; usar esto
como guía para explicar cada diapositiva.

## Repartición sugerida

- Rivera Morales David: portada, introducción y datos.
- Lopez Bernal Yeimi Lizet: objetivo y pipeline.
- Hidalgo Carrillo Amir Gilberto: comparación de métodos y resultados.
- Badager Estrada Aaron Omar: interpretación biológica, limitaciones y cierre.

Si presentan menos personas, se puede dividir en dos bloques: contexto y métodos;
resultados y conclusiones.

## Diapositiva 1. Portada

Buenas tardes. Nuestro proyecto se titula **Identificación de subgrupos clínicos
de COVID-19 mediante expresión génica**. Trabajamos con el dataset público
GSE157103, que contiene datos de RNA-seq de sangre completa. La idea general fue
usar expresión génica para ver si los pacientes con COVID-19 forman grupos que
se parezcan al estado clínico UCI y no UCI.

## Diapositiva 2. Introducción

COVID-19 severo no es un solo estado biológico. Dos pacientes pueden tener la
misma infección, pero respuestas muy distintas: inflamación, respuesta inmune,
daño vascular, comorbilidades y diferentes decisiones clínicas.

En el conjunto local quedaron 126 muestras alineadas con metadatos: 100 muestras
COVID-19 y 26 controles. La etiqueta clínica principal para evaluar el
agrupamiento fue UCI contra no UCI.

## Diapositiva 3. Datos verificados

Después del filtrado quedaron 11,924 genes. Para el análisis de agrupamiento
usamos los 3,000 genes más variables, porque esos genes concentran más señal
entre pacientes.

La expresión se transformó como `log2(TPM + 1)` y después se escaló por gen. Esto
ayuda a que genes con valores muy grandes no dominen todo el análisis.

También verificamos que GSE157103 está publicado en GEO/NCBI y que los archivos
locales de conteos y TPM contienen las 126 muestras alineadas.

## Diapositiva 4. Pipeline

El pipeline tuvo cinco pasos importantes. Primero se filtraron genes con baja
señal. Después se transformó la expresión con `log2(TPM + 1)`. Luego usamos PCA
para resumir miles de genes en componentes principales.

Con esa representación aplicamos los métodos de agrupamiento. La etiqueta UCI/no
UCI no se usó para construir los grupos; se usó al final, solo para medir qué
tanto coincidían los grupos con una variable clínica real.

## Diapositiva 5. Objetivo

El objetivo fue evaluar si un análisis no supervisado de expresión génica puede
identificar grupos moleculares parcialmente relacionados con UCI/no UCI.

La hipótesis fue que la sangre contiene una señal transcriptómica de severidad,
aunque no esperábamos una separación perfecta. En datos clínicos reales es normal
que haya mezcla, porque la severidad depende de muchos factores.

## Diapositiva 6. Métodos comparados

Se probaron tres familias de métodos. K-means se evaluó con 5 valores de `k`.
El agrupamiento jerárquico se evaluó con 15 combinaciones, incluyendo Ward,
average y complete. DBSCAN se evaluó con 44 combinaciones usando una búsqueda de
parámetros basada en distancias k-NN.

Esto es importante porque DBSCAN sí aparece en la comparación final. No fue el
mejor método: su mejor ARI fue 0.060 y marcó bastante ruido. Por eso no se eligió
como resultado principal.

## Diapositiva 7. Mejor resultado

El mejor resultado fue agrupamiento jerárquico con enlace Ward y `k=2`.

El ARI fue 0.263, la silueta fue 0.176, el NMI fue 0.207 y la exactitud binaria
aproximada fue 0.76. Estos valores no significan clasificación perfecta. El ARI
mide qué tanto coinciden los grupos con UCI/no UCI, y la silueta mide qué tan
separados están internamente los grupos.

La lectura correcta es que existe una señal parcial, no una división clínica
perfecta.

## Diapositiva 8. Grupo vs. UCI

Esta tabla muestra la relación entre los grupos encontrados y el estado clínico.
El grupo 0 tiene 40 pacientes UCI y 14 no UCI. En cambio, el grupo 1 tiene 10
pacientes UCI y 36 no UCI.

La observación principal es que el grupo 0 concentra más pacientes UCI, mientras
que el grupo 1 concentra más pacientes no UCI. Hay mezcla, pero la dirección es
clara. Por eso decimos que la señal existe, aunque no funciona como diagnóstico
perfecto.

## Diapositiva 9. Genes marcadores

Después revisamos genes aumentados por grupo. En el grupo 0 aparecen genes como
`IL1R2`, `UPP1`, `OLAH` y `LRG1`, que son compatibles con activación inmune
innata e inflamación.

En el grupo 1 aparece `CD4` y otros genes que sugieren un perfil más relacionado
con células T y regulación celular.

Esto ayuda a interpretar los grupos: no solo son puntos separados en PCA, también
tienen diferencias biológicas que se pueden leer.

## Diapositiva 10. Enriquecimiento funcional

El enriquecimiento funcional fue más claro para el grupo 0. Los términos
principales fueron degranulación de neutrófilos, respuesta inflamatoria,
respuesta de defensa y sistema inmune innato.

Esta parte es importante porque conecta el agrupamiento con biología. El grupo
con más pacientes UCI también muestra procesos esperados en enfermedad severa:
inflamación y respuesta inmune innata.

La interpretación sigue siendo exploratoria. No estamos diciendo que estos genes
diagnostiquen por sí solos a un paciente, sino que ayudan a explicar la señal
molecular observada.

## Diapositiva 11. Limitaciones

El análisis tiene varias limitaciones. Primero, sangre completa mezcla
composición celular y activación transcripcional. Si un grupo tiene más
neutrófilos, parte de la señal puede venir de composición celular.

Segundo, no ajustamos por edad, sexo, ventilación mecánica ni comorbilidades.
Tercero, los marcadores se calcularon de forma exploratoria sobre `log2(TPM +
1)`, no con un modelo formal como DESeq2 o limma-voom.

Aun así, para un análisis exploratorio, los resultados son coherentes y
reproducibles.

## Diapositiva 12. Conclusiones

La conclusión principal es que la expresión génica en sangre sí recupera
parcialmente el estado UCI/no UCI en GSE157103.

El mejor modelo fue jerárquico Ward con `k=2`, con ARI de 0.263 y exactitud
binaria aproximada de 0.76. La parte más fuerte del proyecto no es la
clasificación perfecta, sino la coherencia entre tres cosas: grupos, genes
marcadores y procesos inflamatorios.

Como trabajo futuro, convendría validar los marcadores con DESeq2 o limma-voom y
ajustar por variables clínicas como edad, sexo, ventilación y comorbilidades.

## Diapositiva 13. Bibliografía

Las fuentes principales fueron el estudio de Overmyer et al. en *Cell Systems*,
la ficha oficial de GEO/NCBI para GSE157103, scikit-learn para los métodos de
agrupamiento, UMAP para visualización y g:Profiler para enriquecimiento
funcional.

También se usó Benjamini-Hochberg para controlar falsos descubrimientos en el
análisis exploratorio de genes marcadores.

## Cierre corto

En resumen, el proyecto muestra que los transcriptomas sanguíneos contienen una
señal molecular relacionada con severidad de COVID-19. La señal no sustituye una
evaluación clínica, pero sí permite generar hipótesis biológicas sobre
inflamación, respuesta innata y diferencias entre pacientes UCI y no UCI.

## Posibles preguntas y respuestas

### ¿Por qué no ganó DBSCAN?

DBSCAN depende mucho de densidad y de los parámetros `eps` y `min_samples`. En
este conjunto marcó muchas muestras como ruido o formó grupos con baja
concordancia clínica. Por eso quedó documentado, pero no fue el mejor resultado.

### ¿Por qué usar PCA antes de agrupar?

La matriz original tiene miles de genes. PCA resume la variación principal y
reduce ruido antes de aplicar los métodos de agrupamiento.

### ¿Por qué el ARI no es más alto?

UCI/no UCI es una etiqueta clínica compleja. Depende de biología, pero también de
edad, comorbilidades, tratamiento, momento de muestreo y criterios clínicos. Por
eso una señal parcial es razonable.

### ¿El modelo sirve para diagnosticar severidad?

No. El análisis es exploratorio. Sirve para ver si hay subgrupos moleculares y
para generar hipótesis, no para diagnosticar pacientes de forma individual.

### ¿Qué mejorarían si hubiera más tiempo?

Haría expresión diferencial formal con DESeq2 o limma-voom, ajustaría por
covariables clínicas y evaluaría estabilidad de clusters con remuestreo.
