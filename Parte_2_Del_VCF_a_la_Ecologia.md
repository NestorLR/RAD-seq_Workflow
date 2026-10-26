# Parte II. Del VCF a la ecología: análisis genómico de *Erigeron canadensis*

## Objetivos de aprendizaje

En esta segunda parte del taller utilizaremos los datos genómicos y ambientales publicados para *Erigeron canadensis* para explorar cómo la variación genética puede relacionarse con la distribución geográfica, la historia de invasión y las condiciones ambientales de una especie.

A diferencia de la Parte I, donde procesamos lecturas y obtuvimos
variantes mediante Stacks, aquí partiremos de un VCF publicado. Por
tanto, esta práctica se centra en análisis posteriores al llamado y
filtrado de SNPs.

Al finalizar, podrás:

-   Importar un VCF y comprobar su correspondencia con los metadatos.
-   Explorar la estructura genética mediante PCA.
-   Visualizar e interpretar las proporciones de ancestría
    publicadas.
-   Estimar estadísticas de diversidad genética por población y
    resumirlas por región.
-   Resumir gradientes climáticos mediante una PCA ambiental.
-   Explorar la asociación entre variación genética y predictores
    ambientales mediante dbRDA.
-   Distinguir entre patrones descriptivos, asociaciones estadísticas e
    hipótesis sobre procesos evolutivos.


## 1. Introducción: de los SNPs a las preguntas ecológicas

La genómica de poblaciones permite estudiar cómo se distribuye la
variación genética dentro y entre poblaciones. Sin embargo, miles de
SNPs no constituyen por sí mismos una explicación biológica: para
interpretar los patrones se necesita considerar también la geografía, el
ambiente y la historia demográfica.

*Erigeron canadensis* es una planta anual nativa de Norteamérica que se
ha establecido ampliamente fuera de su rango nativo, especialmente en
Eurasia. Lucas et al. (2026) analizaron poblaciones distribuidas en el
hemisferio norte para estudiar la estructura genética, la diversidad, la
conectividad y la asociación entre genética y ambiente.

Las preguntas que guiarán esta práctica son:

1.  ¿Cómo se distribuye la variación genética entre los individuos?
2.  ¿Qué patrones de ancestría se comparten entre los rangos nativo y no
    nativo?
3.  ¿Cómo varía la diversidad genética entre poblaciones y regiones?
4.  ¿Qué gradientes resumen las condiciones climáticas de los sitios
    muestreados?
5.  ¿Existe asociación estadística entre variación genética y
    predictores ambientales o geográficos?
6.  ¿Qué explicaciones son compatibles con los resultados y cuáles
    requieren análisis adicionales?


## 2. Datos de partida y organización del proyecto

### 2.1. Datos genómicos publicados

El artículo describe un conjunto de 640 individuos de 280 poblaciones
agrupadas en 19 regiones. Después de los filtros, los autores informan
11,501 SNPs distribuidos en 5,080 loci, reteniendo un SNP por locus para
los análisis descritos en el artículo.

El VCF utilizado en esta práctica contiene 5,080 registros de variantes
y 640 muestras. Es importante reportar las dimensiones del archivo que
realmente se está analizando y no confundir el número de registros del
VCF con el número de SNPs mencionado en otras etapas del procesamiento.

### 2.2. Organización de los archivos

Los archivos de entrada utilizados en los scripts se encuentran en
`data/published/`:

``` text
RAD-seq_Workflow/
└── data/
    └── published/
        ├── Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf
        └── Erigeron_canadensis_metadata_students.csv
```

### 2.2. Metadatos de los individuos

Además de la información genómica, utilizaremos algunos metadatos asociados a cada individuo. Estos incluyen información geográfica, clasificación del rango nativo o no nativo y variables ambientales.

Entre las variables disponibles se encuentran:

- Identificador del individuo.
- Rango de distribución (`Range`).
- Región geográfica (`Region_stat`).
- Latitud y longitud.
- Variables bioclimáticas (`Bio01`–`Bio19`).
- Déficit hídrico climático (`climatic_water_deficit`).
- Huella humana (`human_footprint`).
- Componentes de ancestría publicados (`Admixture1`–`Admixture4`).

Revisemos la tabla de metadata en R:

```R
# Importar metadatos
metadata <- read.csv(
  "data/processed/Erigeron_canadensis_metadata_students.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Dimensiones y estructura
dim(metadata)
str(metadata)

# Revisar las primeras filas
head(metadata)

# Número de individuos y duplicados
nrow(metadata)
sum(duplicated(metadata$Ind_name_genind))
sum(is.na(metadata$Ind_name_genind))

# Individuos por rango y región
table(metadata$Range, useNA = "ifany")
table(metadata$Region_stat, useNA = "ifany")

# Resumen de las variables bioclimáticas
summary(metadata[, bio_columns])
```

### 2.4. Inspección del VCF

Ejecuta:

``` bash
Rscript bin/05_vcf_check.R
```

El script comprueba las dimensiones del VCF, los identificadores de las
muestras y los datos faltantes. Antes de cualquier análisis, verifica
que los individuos de la matriz genética y la metadata estén
correctamente emparejados y en el mismo orden.


## 3. Estructura genética: análisis de componentes principales

Ejecuta:

``` bash
Rscript bin/06_genetic_PCA.R
```

El PCA genético resume la variación genómica en ejes que representan
direcciones de variación entre individuos. Los primeros componentes
explican la mayor proporción de la variación representada por la matriz
utilizada, pero no necesariamente toda la estructura genética relevante.

### Qué resultados se generarán y cómo interpretarlos

Examina el porcentaje de varianza explicado por PC1 y PC2 y la
distribución de los individuos coloreados por rango o región. Si los
rangos se superponen en los primeros componentes, eso indica que no
existe una separación evidente en esa proyección; no demuestra ausencia
de estructura genética. Patrones asociados con regiones o linajes pueden
aparecer en componentes posteriores o requerir métodos diferentes.

Lucas et al. (2026) describen estructura genética a escala geográfica y
cuatro grupos de ancestría, pero su interpretación se apoyó en varios
análisis, no únicamente en un PCA. Por ello, no se debe esperar que el
PCA por sí solo reproduzca todas las agrupaciones descritas en el
artículo.

### Preguntas de interpretación

-   ¿Qué porcentaje de la variación explican PC1 y PC2?
-   ¿Se distinguen agrupamientos o existe solapamiento?
-   ¿Cambian los patrones al colorear por rango y por región?
-   ¿Qué limitaciones tiene interpretar únicamente los dos primeros
    componentes?

## 4. Ancestría genética y estructura poblacional

Ejecuta:

``` bash
Rscript bin/07_ADMIXURE_visualization.R
```

Este script visualiza las proporciones de ancestría publicadas para
`K = 4`; no vuelve a ejecutar ADMIXTURE sobre todos los individuos. Cada
componente representa una proporción estimada por el modelo y no debe
interpretarse automáticamente como una especie, población biológica
discreta o linaje con historia evolutiva completamente independiente.

### Qué reportaron los autores

Lucas et al. (2026) identificaron `K = 4` como el valor favorecido por
el error de validación cruzada, aunque señalaron que el gráfico no
presentaba un codo inequívoco. Los cuatro componentes mostraron patrones
geográficos, pero tres aparecieron en ambos rangos. Un cuarto componente
fue poco frecuente en el rango nativo y más extendido en el no nativo.
Los autores interpretaron esta distribución como compatible con una
historia compleja de diferenciación, introducciones múltiples y mezcla
genética, sin considerar que las proporciones de ancestría por sí solas
resuelvan el origen histórico de cada población.

### Qué examinar en nuestros gráficos

Compara las barras individuales y las proporciones promedio por
población o región. Observa si los componentes aparecen en ambos rangos,
si hay poblaciones con ancestría mixta y si la composición varía
geográficamente. El promedio por población facilita la visualización
regional, pero puede ocultar variación entre individuos de una misma
localidad.

### Preguntas de interpretación

-   ¿Qué componentes aparecen en ambos rangos?
-   ¿Qué tan homogénea es la composición dentro de cada región?
-   ¿Qué información adicional se necesitaría para inferir las fuentes
    de introducción?
-   ¿Por qué `K = 4` no equivale necesariamente a cuatro poblaciones
    naturales verdaderas?

## 5. Diversidad genética

Ejecuta:

``` bash
Rscript bin/08_genetic_diversity.R
```

En esta sección se resumen estimadores de diversidad por población y,
posteriormente, por región. La unidad de análisis debe estar explícita:
los valores individuales, poblacionales y regionales no son
intercambiables.

Los indicadores incluidos son heterocigosidad observada (`Ho`),
heterocigosidad esperada (`He`), coeficiente de endogamia (`Fis`) y
riqueza alélica (`Ar`).

### Qué reportaron los autores

Lucas et al. (2026) encontraron valores bajos de diversidad a escala poblacional y valores relativamente altos de `FIS`. No detectaron diferencias significativas entre rangos nativo y no nativo para `Ho`, `He`, `FIS` o riqueza alélica mediante las pruebas t descritas en el artículo. Este resultado muestra que, en sus datos y con esos estimadores, no se detectaron diferencias estadísticas entre rangos; no demuestra que ambos rangos sean idénticos ni descarta efectos demográficos de las introducciones.

### Cómo interpretar nuestros resultados

Compara las medias y la dispersión entre regiones, considerando el
número de poblaciones y muestras que contribuyen a cada resumen. No
atribuyas automáticamente valores bajos de diversidad a un efecto
fundador, ni valores de `FIS` a endogamia sin considerar el sistema
reproductivo, la estructura poblacional y el método de estimación.


### Preguntas de interpretación

-   ¿Qué mide cada estadística y a qué nivel se calcula?
-   ¿Qué diferencias descriptivas aparecen entre regiones?
-   ¿Las diferencias visuales se han evaluado estadísticamente?
-   ¿Qué factores, además de la introducción, pueden afectar estos
    indicadores?

## 6. Variación ambiental: PCA de variables climáticas

Ejecuta:

``` bash
Rscript bin/09_enviromental_PCA.R
```

### Qué resultados se generarán y cómo interpretarlos

El gráfico de individuos permite comparar las condiciones climáticas asociadas con las muestras, y el gráfico de variables ayuda a identificar qué variables contribuyen a cada gradiente. La dirección de los vectores y sus relaciones deben interpretarse junto con las cargas y los porcentajes explicados.

### Qué reportaron los autores

Lucas et al. (2026) realizaron una PCA con las 19 variables
bioclimáticas y describieron un PC1 que explicaba el 36% de la variación
climática, asociado principalmente con estacionalidad térmica y
continentalidad. El PC2 explicó el 25% y resumió variación relacionada
con humedad y clima de la estación cálida. Retuvieron PC1 para los
análisis posteriores con el fin de representar el gradiente
macroclimático dominante y reducir colinealidad, mientras que el déficit
hídrico se incorporó como predictor adicional.

### Preguntas de interpretación

-   ¿Cuánta variación explican PC1 y PC2?
-   ¿Qué variables contribuyen más a cada eje?
-   ¿Se observa solapamiento climático entre los rangos?
-   ¿Qué muestras se excluyeron por valores faltantes y cómo podría
    afectar esto la representación?

## 7. Relación entre variación genética y ambiente

Ejecuta:

``` bash
Rscript bin/10_genetic_env_RDA.R
```

La dbRDA relaciona una matriz de distancias genéticas con predictores
ambientales y geográficos. En el script adaptado se utilizan déficit
hídrico, latitud, PC1 ambiental y huella humana. Los círculos de
ancestría superpuestos en la figura son una capa de visualización: las
proporciones de `Admixture` no son la respuesta directa del modelo si la
fórmula usa la matriz genética.

### Qué reportaron los autores

Lucas et al. (2026) encontraron que la latitud, el déficit hídrico y los
componentes climáticos se asociaban con la estructura genética; el
vector de huella humana fue relativamente corto en la figura publicada.
Los dos primeros ejes restringidos explicaron el 55.3% y el 32.5% de la
variación restringida, equivalentes al 12.5% y 7.3% de la variación
genética total, respectivamente. Al superponer la ancestría, los autores
observaron patrones compatibles con asociaciones entre algunos
componentes y gradientes geográficos o climáticos. Interpretaron estos
patrones en el contexto de la historia demográfica y la dispersión, sin
atribuir toda la estructura a la selección ambiental.

### Cómo interpretar nuestro resultado

Reporta la significancia global, el `R²` y el `R²` ajustado. Las pruebas
por término dependen del método utilizado: las pruebas secuenciales
dependen del orden de los predictores, mientras que las marginales
evalúan cada término condicionado a los demás. Revisa la colinealidad
entre predictores y la construcción de la matriz genética antes de
interpretar los vectores como asociaciones biológicas independientes.

### Preguntas de interpretación

-   ¿Es significativo el modelo global?
-   ¿Qué proporción de variación explica el modelo, antes y después del
    ajuste?
-   ¿Qué variables siguen asociadas con la respuesta en las pruebas
    marginales?
-   ¿Qué parte del patrón podría estar asociada con geografía, historia
    demográfica o dispersión?
-   ¿Por qué una asociación genética-ambiental no demuestra por sí sola
    adaptación local?

## 8. Integración de resultados

Los análisis describen aspectos complementarios del conjunto de datos:

-   El PCA genético resume las principales direcciones de variación
    genética.
-   Las proporciones de ancestría permiten visualizar componentes
    genéticos compartidos.
-   Las estadísticas de diversidad describen variación dentro de
    poblaciones y su distribución regional.
-   El PCA ambiental resume gradientes climáticos.
-   La dbRDA evalúa la asociación entre la matriz genética y los
    predictores seleccionados.

Lucas et al. (2026) interpretaron sus resultados en conjunto: la
presencia de componentes de ancestría compartidos entre rangos, la
estructura geográfica y las asociaciones ambientales fueron compatibles
con una historia de múltiples introducciones, mezcla genética y
diferenciación espacial. También señalaron que la relación entre
genética y ambiente difería entre rangos y que la dispersión mediada por
humanos podría debilitar los patrones ambientales en las regiones no
nativas. Esas conclusiones proceden de varios análisis del artículo,
incluidos análisis de distancia genética y flujo génico que no
reproducimos en esta práctica.

Por ello, nuestros ejercicios permiten evaluar algunos patrones, pero no
reproducen toda la evidencia utilizada por los autores para reconstruir
la historia de invasión.

### Preguntas para la discusión

1.  ¿Qué patrones aparecen de manera consistente en más de un análisis?
2.  ¿Qué resultados difieren entre el rango nativo y el no nativo?
3.  ¿Qué inferencias se sostienen directamente con los análisis
    realizados?
4.  ¿Qué procesos históricos son hipótesis compatibles, pero no
    demostradas, por nuestros gráficos?
5.  ¿Qué análisis adicionales serían necesarios para evaluar aislamiento
    por distancia, aislamiento por ambiente o rutas de dispersión?


## Referencia principal

Lucas, M. S., Rosche, C., Hensen, I., Michalski, S. G., Nagy, D. U.,
Gamba, D., Onstein, R. E., Abidkulova, K., Al-Gharaibeh, M., et
al. (2026). *Population genomics of a cosmopolitan weed provides
insights into its local adaptation and recent demographic history*.
Molecular Ecology. https://doi.org/10.1111/mec.70368
