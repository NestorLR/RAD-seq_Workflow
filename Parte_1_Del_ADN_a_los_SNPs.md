# Parte 1. Del DNA a los SNPs

## De una muestra de DNA a una matriz de variantes genéticas

En esta primera parte de la práctica seguiremos el flujo general de un experimento de **ddRAD-seq**, desde las lecturas de secuenciación hasta la obtención de un archivo **VCF** con variantes genéticas.

El objetivo no es solamente aprender a ejecutar programas bioinformáticos, sino comprender **qué información entra y qué información sale de cada etapa**.

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
Identificación de loci
 ↓
Llamado de SNPs
 ↓
Filtrado
 ↓
VCF

Utilizaremos los datos generados en el artículo de [Lucas et al., 2026](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.70368) para replicar parcialmente los principales resultados. Este estudio analizó las poblaciones de **Erigeron canadensis** (Asteraceae) nativas e invasoras usando genómica de poblaciones y genómica del paisaje.  

El diseño de muestreo consisitió en 640 indivudos de 280 poblaciones de 19 regiones (8 nativas y 11 no-nativas) distribuídas a lo largo del Hemisferio Norte. 

## 1.1. Genómica de poblaciones y genómica del paisaje


## 1.2. ¿Qué es RAD-seq?
RAD-seq (Restriction site Associated DNA sequencing) es una estrategia de secuenciación que permite obtener información de muchas regiones distribuidas a lo largo del genoma.
La idea general es:
1. El DNA se corta utilizando enzimas de restricción.
2. Se seleccionan fragmentos asociados a los sitios de corte.
3. Se añaden adaptadores.
4. Se incorporan identificadores para distinguir individuos.
5. Los fragmentos se secuencian.
6. Las lecturas se procesan para identificar loci y variantes.
Esto permite obtener miles de marcadores distribuidos por el genoma sin necesidad de secuenciar el genoma completo.

## 1.3. ¿Qué es ddRAD-seq?
En esta práctica trabajaremos específicamente con ddRAD-seq (double-digest RAD-seq).
A diferencia de otros enfoques RAD-seq, ddRAD utiliza una combinación de dos enzimas de restricción para generar los fragmentos.

*Pregunta*
¿Por qué utilizar solamente una fracción del genoma puede ser útil para estudios de genética de poblaciones?

## 1.4. De la biblioteca de ddRAD a las lecturas FASTQ
Después de construir la biblioteca ddRAD, los fragmentos son secuenciados.
En nuestro caso utilizaremos secuenciación Illumina.

## 1.5. Obteneción de los datos
Los datos crudos de la secuenciación masiva están disponibles en el NCBI  bajo el número de acceso= PRJEB91610. 

```bash
cd bin/
```

```bash
bin/01-Get-data.sh
```

Inspeccciona los scripts con 'nano'

```bash
nano 01-Get-data.sh
```

Para esta parte, solo descargaremos los archivos correspondientes a dos de los 6individuos con el fin de  utilizaremos un subconjunto de los datos originales para que el procesamiento pueda realizarse durante una sesión de clase.


El secuenciador produce archivos FASTQ.
Un archivo FASTQ contiene cuatro líneas por lectura:

```
@READ_ID
ACGTGCTAGCTAGCTAGCTA
+
FFFFFFFFFFFFFFFFFFFF
```
Estas cuatro líneas representan:
1. Identificador de la lectura.
2. Secuencia.
3. Separador.
4. Calidad de cada base.


Inspeccionemos un archivo 

```bash
zcat data/raw/sample_R1.fastq.gz | head -n 20
```

## 1.4. Control de calidad de los datos
Antes de procesar las lecturas debemos evaluar su calidad.
Para ello utilizaremos FastQC y MultiQC.
FastQC permite revisar diferentes características de las lecturas, entre ellas:
- calidad por posición;
- contenido de GC;
- secuencias sobre-representadas;
- presencia de adaptadores;
- distribución de la longitud de las lecturas.

MultiQC hace un reporte unificado de todas las muestras. 

```bash
bin/02-Quality.Control.sh
```
Revisa el reporte de calidad. Si todo está bien,procedemos al ensamble y obtención de variantes. 

## 1.5. Ensamblado de loci y llamado de variantes
Para obtener nuestros datos genéticos por individuos, usaremos el programa Stacks.

### 1.5.1 Paso 1

```bash
bin/03-a-Stacks_process_radtags.sh
```

### 1.5.2 Paso 2

Los archivos procesados en el paso anterior, mostraron una incompatibilidad en el paso siguiente, el cual se resolvió normalizando los encabezados de cada read con el comando awk

```bash
mkdir -p data/stacks/process_radtags_fixed

for sample in ERR15191538 ERR15191667; do
    for mate in 1 2; do
        zcat "data/stacks/process_radtags/${sample}.${mate}.fq.gz" |
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
Posteriormente, procedimos con el ensamble y llamdo de variantes

```bash
bin/03-b-Stacks_denovo-assembly.sh 
```

### 1.5.3 Paso 3

Generamos el archivo VCF
```bash
bin/03-c-Stacks_populations.sh
```