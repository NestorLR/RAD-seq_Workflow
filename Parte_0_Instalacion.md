# Parte 0. Instalación y configuración del ambiente

En esta práctica trabajaremos con datos de **ddRAD-seq** para explorar cómo los datos genómicos pueden utilizarse para responder preguntas ecológicas y evolutivas.

Antes de comenzar con los análisis es necesario instalar las herramientas que utilizaremos durante la práctica.

---

## 0.1. Requisitos

Necesitarás:

- Linux, macOS o Windows + WSL2
- Conda o Mamba
- Conexión a internet
- Al menos ~10 GB de espacio libre
- Terminal

Durante la práctica utilizaremos principalmente:

| Programa | Función |
|---|---|
| FastQC | Evaluar la calidad de las lecturas FASTQ |
| MultiQC | Resumir los resultados de control de calidad |
| Stacks | Procesamiento de datos RAD-seq y generación de SNPs |
| R | Análisis estadísticos y visualización |
| RStudio (opcional) | Interfaz gráfica para trabajar en R |

---

## 0.2. Verificar que Conda está instalado

Abre una terminal y ejecuta:

```bash
conda --version
```

## 0.3. Crear y activar un nuevo ambiente de conda

```bash
conda create -n radseq -c conda-forge python=3.11
```

```bash
conda activate radseq
```

## 0.4 Instalar las herramientas bioinformáticas

```bash
conda install -c conda-forge -c bioconda \
    sra-tools \
    fastqc \
    multiqc \
    stacks \
    r-base
```
y

```bash
conda install -c conda-forge -c bioconda stacks=2.68
```

## 0.5 Instalar paqueterías de R

```bash
conda install -c conda-forge \
    r-vcfR \
    r-adegenet \
    r-hierfstat \
    r-vegan \
    r-ggplot2 \
    r-dplyr \
    r-tidyr \
    r-readr \
    r-stringr \
    r-scatterpie \
    r-maps \
    r-ape
```

Comprobar que las paqueterías funcionan. Entra a R

```bash
R

```

```R
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(vcfR)
library(adegenet)
library(hierfstat)
library(vegan)

```
cuando termine la instalación, salir de R

```R
q()

```

**NOTA:** 
En mi computadora, no pude cargar adegent la primera vez debido a una discordancia con la versión de una paquetería necesaria. Tuve que instalar la versión que el sistema requería:
```bash
conda install -c conda-forge "r-stringi=1.7.12"
```

## Fin. 