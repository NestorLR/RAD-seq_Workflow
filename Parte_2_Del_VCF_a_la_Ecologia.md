# Parte II. Del VCF a la ecología: análisis genómico de *Erigeron canadensis*

## Objetivos de aprendizaje

En esta segunda parte del taller utilizaremos los datos genómicos y ambientales publicados para *Erigeron canadensis* para explorar cómo la variación genética puede relacionarse con la distribución geográfica, la historia de invasión y las condiciones ambientales de una especie.

A diferencia de la Parte I, donde aprendimos a procesar lecturas de secuenciación y obtener variantes genéticas mediante un flujo bioinformático con Stacks, aquí partiremos de la matriz genómica completa proporcionada por los autores del estudio.

Al finalizar esta práctica, serás capaz de:

- Importar una matriz de SNPs en R y relacionarla con los metadatos de los individuos.
- Explorar la estructura genética mediante análisis de componentes principales (PCA).
- Interpretar patrones de ancestría genética y visualizar resultados de ADMIXTURE.
- Calcular e interpretar estadísticas de diversidad genética.
- Resumir la variación ambiental mediante una PCA de variables climáticas.
- Explorar la asociación entre variación genética y condiciones ambientales mediante análisis multivariados.
- Interpretar los resultados desde una perspectiva de ecología evolutiva y biología de invasiones.

## 1. Introducción: de los SNPs a las preguntas ecológicas

La genómica de poblaciones permite estudiar la distribución de la variación genética dentro y entre poblaciones, así como explorar los procesos históricos y ecológicos que contribuyen a generar esos patrones.

Sin embargo, obtener miles de SNPs no constituye por sí mismo una explicación biológica. Para interpretar los patrones genéticos es necesario relacionarlos con información adicional, como la ubicación geográfica de los individuos y las condiciones ambientales.

En esta práctica utilizaremos el estudio de *Erigeron canadensis*, una especie ampliamente distribuida cuyo amplitud nativa se encuentra en Norteamérica y que también se ha establecido fuera de esta región, a Europa.

Trabajaremos con los datos genómicos y los metadatos publicados por los autores para investigar la variación genética entre individuos y explorar su relación con los gradientes geográficos y ambientales.

### Preguntas biológicas

A lo largo de la práctica abordaremos las siguientes preguntas:

1. ¿Cómo se distribuye la variación genética entre los individuos estudiados?
2. ¿Existen patrones de estructura genética o ancestría compartidos entre las regiones nativas y no nativas?
3. ¿Cómo podemos describir y comparar la diversidad genética de los grupos estudiados?
4. ¿Qué gradientes ambientales resumen mejor las condiciones climáticas asociadas con los individuos?
5. ¿Existe una asociación estadística entre la variación genética y las variables ambientales?
6. ¿Qué procesos históricos y ecológicos podrían contribuir a explicar los patrones observados?

**OJO** Estas preguntas se abordarán mediante análisis exploratorios y multivariados. Pero, la existencia de una asociación estadística no demuestra por sí sola que el ambiente haya causado diferenciación genética o que exista adaptación local.

## 2. Datos de partida y organización del proyecto

### 2.1. Matriz genómica publicada

Para esta práctica utilizaremos la matriz genómica completa publicada por los autores, en lugar del pequeño conjunto de dos individuos procesado en la Parte I.

El conjunto de datos final descrito en el artículo comprende 640 individuos y 11,501 SNPs, con un SNP por locus después de los filtros aplicados por los autores.


### 2.3. Organización de los archivos

Los archivos de entrada se encuentran en data/published:

```text
RAD-seq_Workflow/
├── data/
│   ├── published/
│   │   ├── genotypes.vcf
│   │   └── Erigeron_canadensis_metadata_cleaned.csv
│   │   └── Erigeron_canadensis_metadata_students.csv


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


### 2.4. Datos genómicos

Importemos la matriz de SNPs y comprobaremos que los identificadores de los individuos coincidan con la matriz de metadatos.

Revisemos la matriz en R:

También revisaremos:

- El número de individuos y variantes.
- La presencia de identificadores duplicados.
- Los datos faltantes.


```R
Rscript bin/05_vcf_check.R
```

## 3. Estructura genética: análisis de componentes principales

En esta sección utilizaremos una PCA genética para explorar cómo se distribuye la variación genómica entre los individuos.

El PCA resume la variación multivariada en un conjunto reducido de ejes que permiten visualizar patrones de similitud y diferenciación genética.

### Actividades

1. Importar la matriz de SNPs.
2. Revisar y preparar los genotipos.
3. Examinar el tratamiento de los datos faltantes.
4. Calcular los componentes principales.
5. Representar los individuos en los primeros componentes.
6. Colorear los individuos según su rango de distribución.
7. Explorar la relación entre los patrones genéticos y la procedencia geográfica.

### Preguntas de interpretación

- ¿Qué proporción de la variación explica cada componente?
- ¿Se observan agrupamientos genéticos definidos?
- ¿Los individuos del rango nativo y no nativo se separan claramente?
- ¿Qué otras explicaciones podrían generar los patrones observados?


Corre el script 

```R
Rscript bin/06_genetic_PCA.R
```

## 4. Ancestría genética y estructura poblacional

La estructura genética puede explorarse mediante métodos que estiman proporciones de ancestría bajo distintos modelos. En esta sección distinguiremos entre ejecutar un análisis y visualizar resultados publicados.

### Actividades

1. Revisar los fundamentos del análisis de ancestría.
2. Comprender el significado del parámetro K.
3. Ejecutar un ejercicio reducido con ADMIXTURE, cuando el tiempo y los recursos computacionales lo permitan.
4. Visualizar las proporciones de ancestría publicadas para el conjunto completo.
5. Comparar los patrones de ancestría entre regiones nativas y no nativas.


Para visualizar los resultados de ADMIXURE, corre el script 

```R
Rscript bin/07_ADMIXURE_visualization.R
```



### Preguntas de interpretación

- ¿Los individuos presentan componentes de ancestría compartidos?
- ¿Los patrones de ancestría coinciden con los grupos geográficos?
- ¿Qué limitaciones tiene interpretar un valor de K como el número real de poblaciones?
- ¿Qué información adicional sería necesaria para explicar los patrones observados?


## 5. Diversidad genética

En esta sección exploraremos cómo se cuantifica la variación genética y cómo pueden compararse los grupos estudiados.

### Conceptos principales

- Heterocigosidad observada.
- Heterocigosidad esperada.
- Diversidad nucleotídica.
- Coeficiente de endogamia.

### Actividades

1. Revisar las estadísticas de diversidad disponibles en los datos publicados.
2. Identificar la unidad de análisis de cada estadística.
3. Organizar los resultados en tablas.
4. Comparar descriptivamente los grupos.


Es importante distinguir las estadísticas calculadas por individuo de aquellas estimadas por población o grupo. No todas pueden calcularse directamente a partir de una matriz individual sin definir previamente cómo se agruparán las muestras.

Para obtener los estimadores de diversidad genética, corre el script: 

```R
Rscript bin/08_genetic_diversity.R
```

### Preguntas de interpretación

- ¿Qué mide cada estadística de diversidad?
- ¿Qué diferencias se observan entre los grupos?
- ¿Qué información aporta la diversidad genética a los estudios de invasiones biológicas?



## 6. Variación ambiental: PCA de variables climáticas

Para estudiar la relación entre genética y ambiente, primero resumiremos la variación ambiental asociada con los individuos.

Las variables climáticas suelen estar correlacionadas entre sí. Una PCA ambiental permite identificar gradientes que resumen parte de esa variación y facilita la exploración de patrones ambientales sin incluir simultáneamente todas las variables originales.

### Actividades

1. Seleccionar las variables ambientales disponibles.
2. Examinar datos faltantes y distribuciones.
3. Evaluar correlaciones entre variables.
4. Estandarizar las variables cuando corresponda.
5. Calcular el PCA ambiental.
6. Examinar los porcentajes de varianza explicados.
7. Interpretar las cargas de las variables y los gradientes ambientales.
8. Representar los individuos en los principales componentes.

La PCA ambiental se calculará a partir de las variables ambientales originales. Los campos `PC1` y `PC2` de los metadatos que correspondan a la PCA genética no deben confundirse con los componentes ambientales que calcularemos aquí.

Corre el script: 

```R
Rscript bin/09_enviromental_PCA.R
```


### Preguntas de interpretación

- ¿Qué gradientes ambientales explican una mayor proporción de la variación?
- ¿Qué variables contribuyen más a cada componente?
- ¿Se observa diferenciación ambiental entre los rangos nativo y no nativo?
- ¿Qué limitaciones tiene resumir varias variables en unos pocos componentes?


## 7. Relación entre variación genética y ambiente

En esta sección exploraremos si los patrones genéticos están asociados con variables ambientales y geográficas.

Utilizaremos un análisis multivariado de redundancia basado en distancias (dbRDA).

### Actividades

1. Definir la matriz de respuesta genética.
2. Seleccionar las variables ambientales y geográficas.
3. Construir el modelo multivariado.
4. Evaluar la varianza explicada.
5. Examinar la significancia mediante pruebas de permutación.
6. Visualizar los resultados e interpretar los gradientes asociados con la variación genética.


### Preguntas de interpretación

- ¿Existe una asociación estadística entre la variación genética y los predictores ambientales?
- ¿Qué variables presentan mayor asociación con la ordenación?
- ¿Qué otros procesos podrían explicar los patrones observados?

Corre el script: 

```R
Rscript bin/10_genetic_env_RDA.R
```


## 8. Integración de resultados

En esta sección reuniremos los resultados de los análisis anteriores.

La interpretación no dependerá de una sola figura, sino de la comparación de distintos tipos de evidencia:

- La PCA genética describe los principales patrones de variación.
- El análisis de ancestría permite explorar componentes genéticos compartidos.
- Las estadísticas de diversidad resumen aspectos de la variación genética.
- La PCA ambiental identifica gradientes de las condiciones climáticas.
- La dbRDA explora la asociación entre la variación genética y los predictores seleccionados.

### Preguntas para la discusión

1. ¿Qué patrones genéticos se observan entre los individuos?
2. ¿Cómo se relacionan estos patrones con el rango nativo y no nativo?
3. ¿Qué diferencias ambientales se identifican?
4. ¿Existe evidencia de asociación entre genética y ambiente?
5. ¿Qué resultados son consistentes con el estudio publicado?
6. ¿Qué conclusiones no pueden obtenerse únicamente con estos análisis?

Se distinguirán explícitamente los resultados reproducidos a partir de los datos publicados, los análisis nuevos realizados durante la práctica y las hipótesis que requieren evidencia adicional.

## 9. Conclusiones

Al finalizar, cada estudiante deberá poder explicar cómo una matriz de SNPs puede utilizarse para investigar preguntas de ecología evolutiva.

Las conclusiones deberán resumir:

- Los principales patrones de variación genética.
- Las diferencias de diversidad observadas, cuando corresponda.
- Los gradientes ambientales identificados.
- La evidencia de asociación entre genética y ambiente.
- Las limitaciones de los datos y de los métodos empleados.

No se interpretará una asociación estadística como prueba definitiva de adaptación local, ni se inferirán procesos demográficos específicos a partir de una sola ordenación.
