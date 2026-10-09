# Parte 1. Del DNA a los SNPs

## De una muestra de DNA a una matriz de variantes genéticas

En esta primera parte de la práctica seguiremos el flujo general de un experimento de **ddRAD-seq**, desde las lecturas de secuenciación hasta la obtención de un archivo **VCF** (*Variant Call Format*) que contiene información sobre variantes genéticas.

El objetivo no es solamente aprender a ejecutar programas bioinformáticos, sino comprender **qué información entra, qué proceso se realiza y qué información sale de cada etapa**.

El flujo general será:

```text
DNA
 ↓
Digestión con enzimas de restricción
 ↓
Selección de fragmentos
 ↓
Ligación de adaptadores y códigos de barras
 ↓
Secuenciación Illumina
 ↓
FASTQ
 ↓
Control de calidad
 ↓
Procesamiento de ddRAD-seq
 ↓
Identificación y ensamblaje de loci
 ↓
Llamado de SNPs
 ↓
Filtrado de variantes
 ↓
VCF
```

Utilizaremos los datos generados en el artículo de Lucas et al. (2026), *Population genomics of a cosmopolitan weed provides insights into its local adaptation and recent demographic history*, para reproducir parcialmente los principales resultados.

El estudio analizó individuos de **_Erigeron canadensis_** (Asteraceae), una especie con un amplio rango de distribución nativo en Norteamérica y poblaciones establecidas fuera de este rango. Los autores emplearon herramientas de genómica de poblaciones y genómica del paisaje para investigar sus patrones de variación genética y su historia demográfica.

El diseño de muestreo incluyó 640 individuos de 280 poblaciones distribuidas en 19 regiones: ocho nativas y once no nativas.

En esta primera parte utilizaremos únicamente dos individuos para que el procesamiento pueda realizarse durante una sesión de clase. Este subconjunto servirá para aprender el flujo computacional, pero **no será suficiente para reproducir las conclusiones poblacionales del artículo**.

Artículo: [Lucas et al. (2026)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.70368)

---

## 1.1. Genómica de poblaciones y genómica del paisaje

La **genómica de poblaciones** estudia cómo se distribuye la variación genética dentro y entre poblaciones. Para ello utiliza numerosos marcadores distribuidos a lo largo del genoma, lo que permite investigar procesos como la diferenciación genética, el flujo génico, la historia demográfica y la estructura poblacional.

La **genómica del paisaje** integra la información genética con variables geográficas y ambientales. Su objetivo es explorar cómo factores como la distancia geográfica, el clima, la topografía o la fragmentación del hábitat se relacionan con la distribución de la variación genética.

Por ejemplo, podemos plantear las siguientes preguntas:

- ¿Los individuos de regiones geográficamente cercanas son genéticamente más similares?
- ¿Existen diferencias genéticas entre las regiones nativas y no nativas?
- ¿La variación genética está asociada con gradientes ambientales?
- ¿Qué procesos históricos y ecológicos podrían explicar los patrones observados?

Para abordar estas preguntas necesitamos obtener marcadores genéticos de los individuos y relacionarlos posteriormente con información geográfica y ambiental.

**Idea clave:** los datos genéticos por sí solos no explican todos los procesos ecológicos. Para investigar su relación con el ambiente debemos integrar los genotipos con metadatos de cada individuo.

## 1.2. ¿Qué es RAD-seq?

RAD-seq (*Restriction site-Associated DNA sequencing*) es una estrategia de secuenciación que permite estudiar una fracción del genoma utilizando enzimas de restricción.

La idea general es:

1. El DNA se corta utilizando enzimas que reconocen secuencias específicas.
2. Se seleccionan determinados fragmentos de DNA.
3. Se añaden adaptadores que permiten amplificar y secuenciar los fragmentos.
4. Se incorporan identificadores moleculares, como códigos de barras, para distinguir las muestras cuando corresponda.
5. Los fragmentos se secuencian.
6. Las lecturas se procesan para identificar regiones homólogas del genoma y detectar variantes genéticas.

Esto permite obtener miles de marcadores distribuidos por el genoma sin necesidad de secuenciar el genoma completo de cada individuo.

Una de sus principales ventajas es que permite comparar muchos individuos a un costo generalmente menor que el de la secuenciación genómica completa.

Sin embargo, RAD-seq también tiene limitaciones. La información se obtiene únicamente de las regiones recuperadas mediante el protocolo, y la ausencia de una lectura en un individuo no significa necesariamente que la región no exista en su genoma. Puede deberse, por ejemplo, a diferencias en los sitios de restricción, a la calidad de la secuenciación o a una cobertura insuficiente.

## 1.3. ¿Qué es ddRAD-seq?

En esta práctica trabajaremos específicamente con **ddRAD-seq** (*double-digest RAD-seq*).

A diferencia de otros enfoques RAD-seq que utilizan una sola enzima para definir los sitios de restricción, ddRAD-seq utiliza una combinación de dos enzimas. Esto permite definir los extremos de los fragmentos recuperados y seleccionar una fracción del genoma mediante criterios como el tamaño de los fragmentos.

El procedimiento experimental incluye la digestión del DNA, la selección de fragmentos, la construcción de bibliotecas y la secuenciación.

En nuestro flujo computacional no realizaremos estas etapas de laboratorio. Trabajaremos con las lecturas FASTQ que ya fueron generadas mediante secuenciación.

### Pregunta para reflexionar

¿Por qué puede ser útil estudiar solamente una fracción del genoma en lugar de secuenciar el genoma completo de todos los individuos?

Considera el costo, la cantidad de datos que se obtiene, el número de individuos que podemos analizar y las limitaciones de la información recuperada.

## 1.4. De la biblioteca ddRAD a las lecturas FASTQ

Después de construir la biblioteca ddRAD, los fragmentos seleccionados se secuencian. En el estudio que utilizaremos se empleó secuenciación Illumina.

El secuenciador no produce directamente una tabla de genotipos. Produce **lecturas de secuenciación**, es decir, secuencias de nucleótidos acompañadas de una estimación de la calidad de cada base.

Estas lecturas se almacenan habitualmente en archivos FASTQ.

Como nuestros datos son paired-end, cada fragmento secuenciado puede estar representado por dos lecturas:

- **R1:** lectura obtenida desde uno de los extremos del fragmento.
- **R2:** lectura obtenida desde el extremo opuesto.

Ambas lecturas corresponden al mismo fragmento original, aunque sus secuencias y longitudes pueden diferir.

Esta correspondencia será importante durante el procesamiento bioinformático: si R1 y R2 dejan de estar correctamente emparejadas, algunos programas no podrán reconstruir adecuadamente la información de los fragmentos.

## 1.5. Obtención de los datos

Los datos de secuenciación del estudio están disponibles públicamente en el bioproyeto **PRJEB91610**.

Puedes consultar los registros en [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/view/PRJEB91610).

En esta práctica descargaremos únicamente dos muestras. Esto reduce el tiempo de procesamiento y permite trabajar con los datos durante una sesión de clase.

### 1.5.1. Descargar los datos

Desde la raíz del repositorio, entra al directorio de scripts:

```bash
cd bin/
```

Inspecciona el script de descarga antes de ejecutarlo:

```bash
nano 01-Get-data.sh
```

**¿Por qué inspeccionamos el script?**

Antes de ejecutar un programa, conviene saber qué archivos descargará, dónde los guardará y qué comandos utilizará. Leer los scripts también ayuda a detectar rutas incorrectas o descargas innecesariamente grandes.

Ejecuta el script:

```bash
bash 01-Get-data.sh
```

Regresa a la raíz del repositorio cuando sea necesario:

```bash
cd ..
```

Comprueba que los archivos se hayan descargado:

```bash
ls -lh data/raw/
```

En nuestro conjunto de práctica esperamos encontrar los archivos correspondientes a las dos muestras seleccionadas:

```text
ERR15191538.1.fastq.gz
ERR15191538.2.fastq.gz
ERR15191667.1.fastq.gz
ERR15191667.2.fastq.gz
```

Los sufijos `.1` y `.2` identifican las dos lecturas paired-end de cada muestra.

### 1.5.2. Inspeccionar un archivo FASTQ

Un archivo FASTQ contiene cuatro líneas por lectura:

```text
@READ_ID
ACGTGCTAGCTAGCTAGCTA
+
FFFFFFFFFFFFFFFFFFFF
```

Estas líneas representan:

1. **Identificador:** permite reconocer la lectura.
2. **Secuencia:** las bases de DNA observadas.
3. **Separador:** delimita la secuencia y la información de calidad.
4. **Calidad:** valores codificados que representan la confianza del secuenciador en cada base.


Inspecciona las primeras lecturas de una muestra:

```bash
zcat data/raw/ERR15191538.1.fastq.gz | head -n 20
```

El comando `zcat` permite leer archivos comprimidos con gzip sin descomprimirlos permanentemente. `head -n 20` muestra las primeras 20 líneas, equivalentes a cinco registros FASTQ completos.

Para examinar la lectura correspondiente del otro extremo, cambia `.1` por `.2` en la ruta.

---

## 1.6. Control de calidad de los datos

Antes de procesar las lecturas debemos evaluar su calidad.

Para ello utilizaremos **FastQC** y **MultiQC**.

FastQC genera un informe por archivo que resume distintas características de las lecturas, entre ellas:

- Calidad por posición.
- Distribución de la calidad.
- Contenido de GC.
- Distribución de la longitud de las lecturas.
- Secuencias sobre-representadas.
- Posible presencia de adaptadores.

MultiQC reúne los informes de múltiples archivos en un único reporte, facilitando la comparación entre muestras y entre lecturas R1 y R2.

### 1.6.1. Ejecutar el control de calidad

Inspecciona el script:

```bash
nano bin/02-Quality.Control.sh
```

Ejecuta el análisis:

```bash
bash bin/02-Quality.Control.sh
```

Después, localiza y abre el reporte HTML de MultiQC generado por el script.

### 1.6.2. ¿Cómo interpretar el reporte?

Los resultados de FastQC deben interpretarse considerando el tipo de biblioteca y el diseño experimental.

Por ejemplo, una advertencia sobre contenido de GC o secuencias sobre-representadas no significa automáticamente que la muestra deba descartarse. En datos de RAD-seq pueden aparecer patrones asociados con la propia naturaleza de las bibliotecas.

Debemos identificar si existen problemas que puedan afectar los análisis posteriores, como una pérdida marcada de calidad en los extremos de las lecturas, presencia de adaptadores o diferencias importantes entre archivos.

### Preguntas de comprobación

1. ¿La calidad cambia a lo largo de la longitud de las lecturas?
2. ¿Las lecturas R1 y R2 presentan patrones similares?
3. ¿Hay indicios de adaptadores u otras secuencias sobre-representadas?


**Resultado esperado:** un reporte de calidad que permita decidir si las lecturas pueden continuar al procesamiento bioinformático y si se requiere alguna medida adicional.

---

## 1.7. Procesamiento de lecturas y ensamblaje de loci con Stacks

En esta práctica utilizaremos **Stacks** v2.68 para procesar las lecturas y obtener una matriz de variantes genéticas.

Stacks es un conjunto de herramientas diseñado para trabajar con datos de secuenciación de sitios de restricción, como RAD-seq y ddRAD-seq.

Nuestro flujo utiliza principalmente tres componentes:

- `process_radtags`: procesa las lecturas y aplica filtros iniciales.
- `denovo_map.pl`: coordina el ensamblaje de loci, la construcción del catálogo, la asociación de loci entre muestras y el genotipado.
- `populations`: aplica filtros poblacionales y permite exportar los resultados, por ejemplo, en formato VCF.

En este ejercicio utilizaremos un ensamblaje *de novo*, es decir, **no partiremos de un genoma de referencia para localizar las lecturas**. Los loci se identifican mediante la similitud entre las lecturas y se comparan entre los individuos.

### 1.7.1. Paso 1. Procesamiento de las lecturas con `process_radtags`

Ejecuta el script:

```bash
bash bin/03-a-Stacks_process_radtags.sh
```

Inspecciona el script antes de ejecutarlo para identificar los archivos de entrada, el directorio de salida y los parámetros utilizados.

`process_radtags` procesa las lecturas de secuenciación y puede realizar operaciones como el demultiplexado, la comprobación de sitios de restricción, el filtrado de lecturas de baja calidad y la eliminación de lecturas que no cumplen los criterios establecidos.

En nuestro conjunto de datos, las muestras ya están separadas por individuo, por lo que este paso solo sirve como un pre-procesamiento de los datos.

El programa también debe conservar la correspondencia entre R1 y R2 cuando trabaja en modo paired-end.

#### ¿Qué archivos obtenemos?

El script genera archivos FASTQ procesados y registros de ejecución. El nombre y la ubicación exactos dependen de las opciones especificadas en el script.

Inspecciona los archivos resultantes:

```bash
find data/stacks -maxdepth 2 -type f
```

### 1.7.2. Paso 2. Normalización de los encabezados FASTQ

Durante esta práctica encontramos un problema de compatibilidad al intentar ensamblar los datos paired-end con Stacks 2.68.

Aunque los archivos R1 y R2 tenían el mismo número de lecturas y sus identificadores coincidían, `tsv2bam` no conseguía encontrar las lecturas pareadas correspondientes.

Al revisar los encabezados, observamos que contenían información adicional después del identificador principal. Para este conjunto de datos, normalizamos los encabezados antes de proceder con el ensamblaje.

El siguiente script crea copias de los FASTQ procesados con encabezados simplificados:

```bash
mkdir -p data/stacks_paired_fixed

for sample in ERR15191538 ERR15191667; do
    for mate in 1 2; do
        zcat "data/stacks/${sample}.${mate}.fq.gz" |
        awk -v mate="$mate" '
            NR % 4 == 1 {
                id = $1
                sub(/^@/, "", id)
                sub(/\/[12]$/, "", id)
                print "@" id "/" mate
                next
            }
            { print }
        ' | gzip -c > "data/stacks_paired_fixed/${sample}.${mate}.fq.gz"
    done
done
```

**Antes de ejecutarlo:** comprueba que las rutas de entrada correspondan con los archivos generados por tu script de `process_radtags`. Si el script guarda los FASTQ en otro directorio, modifica las rutas para que coincidan.

#### ¿Qué hace este código?

- `mkdir -p` crea el directorio de salida si todavía no existe.
- `zcat` lee los archivos comprimidos.
- `awk` identifica las líneas de encabezado, que son la primera línea de cada registro FASTQ.
- `$1` extrae el primer campo del encabezado, delimitado por espacios.
- `sub(/^@/, "", id)` elimina el símbolo inicial `@`.
- `sub(/\/[12]$/, "", id)` elimina un sufijo `/1` o `/2` si ya está presente.
- `print "@" id "/" mate` escribe el identificador seguido del sufijo correspondiente al mate.
- `gzip -c` comprime el resultado sin modificar el archivo original.

El script conserva las líneas de secuencia, separador y calidad. La modificación se limita a los encabezados.

Por ejemplo, un encabezado como:

```text
@ERR15191538.1 A00709:533:HLWY5DSX5:4:1101:12319:1031 length=143/1
```

se transforma en:

```text
@ERR15191538.1/1
```

La pareja R2 conserva el identificador principal y termina en `/2`.

Esta normalización resolvió el problema en nuestro conjunto de práctica. Sin embargo, no debe considerarse un paso obligatorio para todos los datos de Stacks: es una solución técnica aplicada a estos archivos y a la incompatibilidad observada.

### 1.7.3. Paso 3. Ensamblaje de loci y genotipado

Una vez preparados los archivos pareados, ejecutaremos el enamblado:

```bash
bash bin/03-b-Stacks_denovo-assembly.sh
```

En este flujo, `denovo_map.pl` coordina varias etapas:

| Programa | Función |
|---|---|
| `ustacks` | Identifica y agrupa lecturas similares dentro de cada individuo para construir loci. |
| `cstacks` | Construye un catálogo y relaciona loci entre individuos. |
| `sstacks` | Asocia los loci de cada individuo con el catálogo. |
| `tsv2bam` | Organiza la información de los loci y recupera las lecturas paired-end correspondientes. |
| `gstacks` | Ensambla las lecturas cuando corresponde y realiza el genotipado. |
| `populations` | Resume la variación y puede aplicar filtros y exportar archivos de variantes. |

El conjunto de parámetros determina qué tan estrictamente se agrupan las lecturas y se relacionan los loci entre individuos. Por ejemplo, `-M` y `-n` controlan distancias de agrupamiento en diferentes etapas del ensamblaje.

No es necesario memorizar todos los parámetros en esta primera práctica. Lo importante es reconocer que el ensamblaje depende de decisiones bioinformáticas que pueden influir en los loci y SNPs recuperados.

#### ¿Qué debemos revisar en el registro?

Al finalizar, comprueba que el registro de `denovo_map.pl` muestre que las etapas principales terminaron correctamente.

También identifica:

- El número de loci ensamblados por individuo.
- El tamaño final del catálogo.
- La proporción de lecturas pareadas recuperadas.
- El número de loci genotipados.
- El número de loci con contigs pareados ensamblados.


**Idea clave:** ensamblar loci, construir un catálogo y llamar genotipos son operaciones distintas. El ensamblaje identifica regiones comparables; el genotipado estima los alelos observados en cada individuo.

---

## 1.8. Exportación y filtrado de variantes: `populations`

Una vez completado el ensamblaje y el genotipado, utilizaremos `populations` para generar archivos que puedan utilizarse en análisis posteriores.

Stacks puede ejecutar `populations` automáticamente como parte de `denovo_map.pl`. Sin embargo, ejecutar `populations` de forma independiente permite especificar filtros y formatos de salida diferentes sin repetir todo el ensamblaje.

En nuestro repositorio utilizaremos un script independiente:

```bash
bash bin/03-c-Stacks_populations.sh
```

Este script toma como entrada los resultados del ensamblaje y el archivo de asignación de muestras a poblaciones (`popMap.txt`).

El archivo de población relaciona cada individuo con el grupo que se utilizará en los análisis. La asignación debe ser consistente con el diseño de muestreo y con las preguntas biológicas.

### 1.8.1. ¿Qué es un VCF?

VCF (*Variant Call Format*) es un formato de archivo utilizado para almacenar variantes genéticas y los genotipos de los individuos.

Entre sus campos se encuentran:

- `CHROM`: identificador de la secuencia o región.
- `POS`: posición de la variante dentro de esa secuencia.
- `REF`: alelo de referencia.
- `ALT`: alelo o alelos alternativos.
- `QUAL`: calidad asociada con la variante.
- `FILTER`: estado de los filtros aplicados.
- `INFO`: información adicional sobre la variante.
- `FORMAT`: descripción de los campos de genotipo.
- Columnas de muestras: genotipos de los individuos.

En un ensamblaje *de novo*, los identificadores de secuencia y las posiciones del VCF corresponden a las secuencias de referencia construidas por el flujo, no necesariamente a posiciones de un cromosoma de un genoma de referencia externo.

Un VCF puede contener genotipos como `0/0`, `0/1` y `1/1` para individuos diploides. Estos representan, respectivamente, dos copias del alelo de referencia, una copia de cada alelo y dos copias del alelo alternativo, siempre que esos genotipos estén codificados de esa forma en el archivo.

### 1.8.2. Filtrado de SNPs

El programa `populations` permite aplicar criterios para decidir qué loci o variantes se conservarán.

Por ejemplo, `-R` establece la proporción mínima de individuos del conjunto global en los que debe estar presente un locus. Otros parámetros permiten controlar la presencia por población, la frecuencia alélica y la profundidad de genotipado.

Los filtros deben elegirse de acuerdo con el diseño de muestreo y los objetivos del análisis. Un filtro demasiado estricto puede eliminar información útil, mientras que uno demasiado permisivo puede conservar variantes con muchos datos faltantes o estimaciones poco confiables.

En nuestro conjunto de dos individuos, los filtros deben interpretarse como un ejercicio técnico. No debemos utilizar estos datos para obtener conclusiones robustas sobre diversidad o diferenciación entre poblaciones.

### 1.8.3. Comprobar el resultado

Después de ejecutar el script, inspecciona el directorio de salida y confirma que el archivo VCF se haya generado.

```bash
find data/stacks/populations_R70 -maxdepth 1 -type f
```

El VCF generado en esta etapa constituye la salida principal de la Parte I. En la Parte II utilizaremos la matriz genómica completa publicada por los autores para desarrollar los análisis ecológicos y de genómica de poblaciones.

---

## 1.9. Síntesis: del DNA a los SNPs

En esta práctica recorrimos las principales etapas de un flujo de análisis de ddRAD-seq:

1. Descargamos lecturas de secuenciación.
2. Inspeccionamos los archivos FASTQ.
3. Evaluamos la calidad de las lecturas.
4. Procesamos las lecturas con `process_radtags`.
5. Preparamos los encabezados paired-end para resolver una incompatibilidad observada.
6. Ensamblamos loci y realizamos el genotipado con Stacks.
7. Exportamos variantes mediante `populations`.

### FIN