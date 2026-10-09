# README — `Erigeron\_canadensis\_metadata\_cleaned.csv`

## Dataset description

Sample-level metadata, geographic coordinates, climate descriptors, greenhouse and field phenotypes, and population-genetic summary statistics for **640 individuals of *Erigeron canadensis* (L.) Cronquist** (syn. *Conyza canadensis*) sampled across its native (North America) and non-native (Eurasia) ranges. Each row is one individual; each column is described below.

* **Format:** comma-separated values (CSV), UTF-8
* **Rows × columns:** 640 × 49
* **Missing data:** empty cells (`NA` when read into R)



To read the file in R:

```r
df <- read.csv("Erigeron\_canadensis\_metadata\_cleaned.csv",
               stringsAsFactors = FALSE,
               na.strings = c("", "NA"),
               check.names = FALSE)
```

\---

## Column descriptions

### Sample identification and location

|#|Column|Description|
|-|-|-|
|1|`Ind\_name\_genind`|Unique individual identifier as used in the accompanying `genind` (adegenet) R object. The population code embedded in this string identifies the sampling locality.|
|2|`Range`|Native/introduced status of the population. Two levels: `Native` (North America) or `Non-native` (Eurasia).|
|3|`Region\_stat`|Biogeographic region used for grouped analyses (19 levels): South-Western North America, Pacific Coast, Intermountain US, Central US, North US Prairies, Great Lakes, South-Eastern US, East Coast, Western Mediterranean, Western Europe, Central Europe, Baltic Europe, Eastern Mediterranean, Eastern Europe, Armenia-Iran, Middle Asia, Siberia, Far East Asia, China.|
|4|`Latitude`|Latitude of the sampling site, decimal degrees, WGS84.|
|5|`Longitude`|Longitude of the sampling site, decimal degrees, WGS84.|
|6|`Collector`|Name of the person who collected the field sample.|

### Greenhouse phenotype

|#|Column|Description|
|-|-|-|
|7|`Initial\_biomass\_from\_greenhouse`|Above-ground dry biomass of the seedling at the start of the common-garden experiment. Units: **grams (g)**.|

### Climatic water balance (1958–2021 mean at the sampling locality)

|#|Column|Description|
|-|-|-|
|8|`climatic\_water\_deficit`|Mean annual climatic water deficit (reference − actual evapotranspiration). Units: **mm yr⁻¹**.|
|9|`CWD\_factor`|Five-level ordinal binning of annual precipitation: `1 <250 mm`, `2 <500 mm`, `3 <750 mm`, `4 <1000 mm`, `5 >1000 mm`.|
|10|`mean\_an\_\_ref\_\_evapotranspiration\_1958-2021`|Mean annual reference (potential) evapotranspiration, 1958–2021. Units: **mm yr⁻¹**.|
|11|`mean\_an\_\_act\_\_evapotranspiration\_1958-2021`|Mean annual actual evapotranspiration, 1958–2021. Units: **mm yr⁻¹**.|
|12|`mean\_an\_\_percipitation\_1958-2021`|Mean annual precipitation, 1958–2021. Units: **mm yr⁻¹**.|



### WorldClim bioclimatic variables (Bio01–Bio19)

Population-level values extracted at the sampling coordinates from **WorldClim v2.1** at **30 arc-second (\~1 km)** resolution (Fick \& Hijmans 2017, *Int. J. Climatol.* 37: 4302–4315).

|#|Column|Description|Unit|
|-|-|-|-|
|13|`Bio01`|Annual mean temperature × 10|°C × 10|
|14|`Bio02`|Mean diurnal temperature range (mean of monthly max−min) × 10|°C × 10|
|15|`Bio03`|Isothermality (Bio02 / Bio07) × 100|%|
|16|`Bio04`|Temperature seasonality (standard deviation × 100)|°C × 100|
|17|`Bio05`|Maximum temperature of warmest month × 10|°C × 10|
|18|`Bio06`|Minimum temperature of coldest month × 10|°C × 10|
|19|`Bio07`|Temperature annual range (Bio05 − Bio06) × 10|°C × 10|
|20|`Bio08`|Mean temperature of wettest quarter × 10|°C × 10|
|21|`Bio09`|Mean temperature of driest quarter × 10|°C × 10|
|22|`Bio10`|Mean temperature of warmest quarter × 10|°C × 10|
|23|`Bio11`|Mean temperature of coldest quarter × 10|°C × 10|
|24|`Bio12`|Annual precipitation|mm|
|25|`Bio13`|Precipitation of wettest month|mm|
|26|`Bio14`|Precipitation of driest month|mm|
|27|`Bio15`|Precipitation seasonality (coefficient of variation)|%|
|28|`Bio16`|Precipitation of wettest quarter|mm|
|29|`Bio17`|Precipitation of driest quarter|mm|
|30|`Bio18`|Precipitation of warmest quarter|mm|
|31|`Bio19`|Precipitation of coldest quarter|mm|

### Field-measured vegetation at the sampling site

Measured in a vegetation plot centred on the focal *E. canadensis* individual.

|#|Column|Description|Unit|
|-|-|-|-|
|32|`Community\_Productivity`|Above-ground dry biomass of all vegetation in the plot.|g m⁻²|
|33|`Conyza\_biomass`|Above-ground dry biomass of *E. canadensis* in the plot.|g|
|34|`Conyza\_biomass\_rate`|Per-day biomass-production rate of *E. canadensis* (biomass divided by plant age in days).|g day⁻¹|
|35|`Conyza\_height`|Maximum height of *E. canadensis* in the plot.|cm|
|36|`Conyza\_inforenscences`|Number of inflorescences of *E. canadensis* per plot. |count|
|37|`Conyza\_cover`|Percent cover of *E. canadensis* in the plot.|%|

### Remote-sensing primary productivity (MODIS Aqua, MYD17A2HGF v6.1)

|#|Column|Description|Unit|
|-|-|-|-|
|38|`GPP\_MYD17A2HGF`|Gross Primary Productivity, 8-day composite, averaged at the sampling coordinates.|kg C m⁻² (8 day)⁻¹|
|39|`PsnNet\_MYD17A2HGF`|Net photosynthesis (PSN\_NET = GPP − maintenance respiration) from the same product.|kg C m⁻² (8 day)⁻¹|

### Population-genetic summary statistics

Derived from the SNP genotype matrix that accompanies this dataset \[Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf].

|#|Column|Description|
|-|-|-|
|40|`He`|Expected heterozygosity (Nei's gene diversity). Unitless, 0–1.|
|41|`Ho`|Observed heterozygosity. Unitless, 0–1.|
|42|`Fis`|Inbreeding coefficient (Wright's *F*<sub>IS</sub>).|
|43|`PC1`|Score on the first principal component of a PCA over the SNP genotype matrix.|
|44|`PC2`|Score on the second principal component of the same PCA.|
|45|`Admixture1`|Ancestry proportion assigned to genetic cluster 1 by ADMIXTURE at K = 4. Sums to 1 across `Admixture1`–`Admixture4`.|
|46|`Admixture2`|Ancestry proportion of cluster 2 (K = 4).|
|47|`Admixture3`|Ancestry proportion of cluster 3 (K = 4).|
|48|`Admixture4`|Ancestry proportion of cluster 4 (K = 4).|

### Human pressure

|#|Column|Description|
|-|-|-|
|49|`human\_footprint`|Human Footprint Index value at the sampling coordinates. Higher values = more human-modified landscapes.|

\---

## Contact

Marilia Souza Lucas — marilia\_souza2@hotmail.com



## Citation

Please cite the associated journal article: DOI: 10.1111/mec.70368

