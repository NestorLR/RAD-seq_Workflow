# Taller: Procesamiento y análisis de datos de tipo *ddRAD-seq* en un contexto de la genómica de poblaciones y genómica del paisaje


Este repositorio contiene los materiales de un taller teórico-práctico de análisis de datos de secuenciación masiva, enfocado en el procesamiento de datos ddRAD-seq y su aplicación en estudios de genética de poblaciones, genómica del paisaje y ecología.

El taller está organizado en dos partes complementarias: la primera consiste en una introducción al procesamiento bioinformático de lecturas de secuenciación de tipo ddRAD-seq, y la segunda, en un flujo de trabajo de análisis de datos genómicos y  ecológicos.

## Objetivos de aprendizaje

Al finalizar el taller, los participantes podrán:

- Comprender los fundamentos de la técnica ddRAD-seq y los principales pasos de su procesamiento bioinformático.
- Realizar controles de calidad y procesar lecturas de secuenciación.
- Comprender el ensamblaje *de novo* de loci y el genotipado con Stacks.

- Explorar patrones de variación genética mediante análisis multivariados.
- Relacionar la variación genómica con la distribución geográfica y las condiciones ambientales.
- Interpretar los resultados considerando sus alcances y limitaciones.

## Contenido del taller

### Parte I. Del ADN a los SNPs

Introducción al procesamiento bioinformático de datos ddRAD-seq, desde las lecturas crudas hasta la identificación de loci y variantes genéticas.



**[Consultar la guía de la Parte I](https://github.com/NestorLR/RAD-seq_Workflow/blob/main/Parte_1_Del_ADN_a_los_SNPs.md)**

### Parte II. Del VCF a la ecología

Análisis de una matriz genómica publicada de *Erigeron canadensis*, integrando información genética, geográfica y ambiental.



**[Consultar la guía de la Parte II](https://github.com/NestorLR/RAD-seq_Workflow/blob/main/Parte_2_Del_VCF_a_la_Ecologia.md)**

## Estructura del repositorio

- `data/`: datos de entrada.
- `metadata/`: información de las muestras y archivos auxiliares.
- `bin/`: scripts para ejecutar los análisis.
- `results/`: resultados, tablas y figuras generadas.
- `README.md`: descripción general y punto de entrada al taller.

La estructura y los nombres de los archivos pueden consultarse directamente en el repositorio.

## Requisitos

Los programas, paquetes y versiones necesarios para cada práctica se describen en la guía de instalación **[guía de instalación](https://github.com/NestorLR/RAD-seq_Workflow/blob/main/Parte_0_Instalacion.md).**

Se recomienda trabajar en Linux o en un entorno compatible con las herramientas de línea de comandos utilizadas en el taller.

## ¿Cómo comenzar?

1. Clona el repositorio.
2. Consulta la guía de instalación y configura el entorno de trabajo.
3. Sigue la Parte I para comprender el procesamiento de lecturas y el ensamblaje de loci.
4. Continúa con la Parte II para analizar los datos genómicos y ambientales.