# RAD-seq_Workflow
Este repositorio documenta el flujo de trabajo bioinformático utilizando SNPs provenientes de datos de secuenciación masiva de tipo RAD-seq



## Estructura del repositorio
El material del curso estará organizado de la sigueinte manera:

RAD-seq_Workflow/
│
├── README.md
│
├── Parte_0_Instalacion.md
├── Parte_1_Del_DNA_a_los_SNPs.md
├── Parte_2_Del_VCF_a_la_ecologia.md
│
├── data/
│   ├── raw/
│   ├── metadata/
│   └── reference/
│
├── bin/
│   ├── 01_fastqc.sh
│   ├── 02_stacks.sh
│   ├── 03_vcf_inspection.sh
│   ├── 04_PCA.R
│   ├── 05_ADMIXTURE.sh
│   ├── 06_ADMIXTURE_plot.R
│   ├── 07_genetic_diversity.R
│   ├── 08_environmental_PCA.R
│   └── 09_RDA.R
│
├── results/
│   ├── expected/
│   └── student/
│
└── figures/
    ├── expected/
    └── paper/