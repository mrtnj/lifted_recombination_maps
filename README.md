# Lifted recombination maps for farm animals


This repo contains recombination maps of farm animal genomes lifted to different
more contemporary genome assemblies.


## Cattle

> Ma, L., O'Connell, J. R., VanRaden, P. M., Shen, B., Padhi, A., Sun, C., ... & Wiggans, G. R. (2015). Cattle sex-specific recombination and genetic control from a large pedigree analysis. PLoS genetics, 11(11), e1005387. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005387

Map in Data Dryad repository: https://datadryad.org/stash/dataset/doi:10.5061/dryad.q2q84

Reference genome: UMD3.1

Lifting:

* UMD3.1.1 -> ARS-UCD1.2


## Chicken

> Groenen, M. A., Wahlberg, P., Foglio, M., Cheng, H. H., Megens, H. J., Crooijmans, R. P., ... & Andersson, L. (2009). A high-density SNP-based linkage map of the chicken genome reveals sequence features correlated with recombination rate. Genome research, 19(3), 510-519. https://genome.cshlp.org/content/19/3/510.short

Map in: Supplemental Table 1 

Reference genome: Galgal2

Supplemental Table 1 was converted to bed, then lifted with UCSC LiftOver (https://genome-euro.ucsc.edu/cgi-bin/hgLiftOver):

Sequential lifting:

* Galgal2 -> Galgal3
* Galgal3 -> Galgal4
* Galgal4 -> Galgal5
* Galgal5 -> Galgal6

