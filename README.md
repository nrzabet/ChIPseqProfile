ChIPseqProfile
==============

R scripts for the http://arxiv.org/abs/1404.5544 paper

Systems requirements

To run the scripts, you need to have access to a machine with an installation of R. To combine multiple ChIP-seq profile plots, the script uses pdflatex command, thus requiring an installation of latex. If the system does not have all required packages (BSgenome.Dmelanogaster.UCSC.dm3 from Bioconductor and xtable) you also need internet connection. Finally, you will need a machine with at least 40 GB of memory and 30 GB of HDD storage space available.
Run the scripts

1. move to the directory. For Linux or MacOS type in the terminal

$ cd ChIPseqProfile

2. start R and then type

$ R

> source("PerformAnalysis.R")

3. (optional) to combine multiple ChIP-seq profile plots open the terminal and call

$ chmod +x GenerateProfilesPanels.sh

$ ./GenerateProfilesPanels.sh

Note that this only works on Linux or MacOS. For Windows systems, please go to the subdirectories listed bellow in section 5 and call pdflatex for the .tex file in each subdirectory
The output

1. The heatmaps can be found in subdirectory

img/heatmap

2. The latex tables with the optimal set of parameters can be found in subdirectory

tables

3. The ChIP-seq profile estimates can be found in subdirectory

img/profile

4. The statistics of the genome wide analysis can be found in subdirectory

img/statistics

5. (optional) The plots with multiple profile files can be found in subdirectories

a). the five ChIP-seq profiles at eve locus using the optimal set of parameters

latex/eve

b). the Bicoid ChIP-seq profiles at all 21 loci using the optimal set of parameters

latex/BCDAccessibleRegions

c). the Caudal ChIP-seq profiles at all 21 loci using the optimal set of parameters

latex/CADAccessibleRegions

d). the Giant ChIP-seq profiles at all 21 loci using the optimal set of parameters

latex/GTAccessibleRegions

e). the Hunchback ChIP-seq profiles at all 21 loci using the optimal set of parameters

latex/HBAccessibleRegions

f). the Kruppel ChIP-seq profiles at all 21 loci using the optimal set of parameters

latex/KRAccessibleRegions

g). the Hunchback ChIP-seq profiles at all 21 loci using the optimal set of parameters and JASPAR PWM motif

latex/HBAccessibleRegionsJASPAR

h). the Kruppel ChIP-seq profiles at all 21 loci using the optimal set of parameters and JASPAR PWM motif

latex/KRAccessibleRegionsJASPAR

i). the Bicoid profile at eve locus using different abundances

latex/BCDEveStripeEnhancer
