# Lifted recombination maps for farm animals


This repo contains recombination maps of farm animal genomes lifted to different
more contemporary genome assemblies.


## Cattle

> Brekke, C., Johnston, S. E., Gjuvsland, A. B., & Berg, P. (2023). Variation and genetic control of individual recombination rates in Norwegian Red dairy cattle. Journal of Dairy Science, 106(2), 1130-1141. https://doi.org/10.3168/jds.2022-22368

Map in Figshare repository https://figshare.com/articles/dataset/Autosomal_linkage_map_NRF/20976067/1.

Reference genome: ARS-UCD1.2

Lifting by BLAT alignment:

* ARS-UCD1.2 -> ARS-UCD1.3



## Chicken

> Groenen, M. A., Wahlberg, P., Foglio, M., Cheng, H. H., Megens, H. J., Crooijmans, R. P., ... & Andersson, L. (2009). A high-density SNP-based linkage map of the chicken genome reveals sequence features correlated with recombination rate. Genome research, 19(3), 510-519. https://genome.cshlp.org/content/19/3/510.short

Map in: Supplemental Table 1, combined sex-averaged map

Reference genome: Galgal2

Supplemental Table 1 was converted to bed, then lifted with UCSC LiftOver (https://genome-euro.ucsc.edu/cgi-bin/hgLiftOver):

Sequential lifting:

* Galgal2 -> Galgal3
* Galgal3 -> Galgal4


> Elferink, M. G., van As, P., Veenendaal, T., Crooijmans, R. P., & Groenen, M. A. (2010). Regional differences in recombination hotspots between two chicken populations. BMC genetics, 11(1), 1-10. https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-11

Map in: Additional file 2

Marker positions in: Additional file 1

Reference genome: Galgal3 (WASHU2)

Map was converted to bed, then lifted with UCSC LiftOver
(https://genome-euro.ucsc.edu/cgi-bin/hgLiftOver):

* Galgal3 -> Galgal4 and Galgal5
* Galgal5 -> GRCg6a