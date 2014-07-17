

#1. install required packacges 
source("InstallPackages.R")

#2. load the libraries 
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm3));
suppressPackageStartupMessages(library(xtable));

source('functions/GenomicGeneralFunctions.R');
source('functions/ComputeOccupancyFunctions.R');
source('functions/GenomicOutputFunctions.R');



#3. process the accessibility data
profile=NULL;
bindingEnergy=NULL;
occupancy=NULL;
DNAAccessibility=NULL;
print("Process accessibility data")
source("ProcessAccessibility.R")

#4. perform the analysis without DNA accessibility
profile=NULL;
bindingEnergy=NULL;
occupancy=NULL;
DNAAccessibility=NULL;
print("Perform the analysis without DNA accessibility")
source("ScriptNoAccessibility.R")


#5. perform the analysis with binary DNA accessibility 
profile=NULL;
bindingEnergy=NULL;
occupancy=NULL;
DNAAccessibility=NULL;
print("Perform the analysis with binary DNA accessibility")
source("ScriptAccessibilityRegions.R")

#6. perform the analysis with binary DNA accessibility and include weak binding sites
profile=NULL;
bindingEnergy=NULL;
occupancy=NULL;
DNAAccessibility=NULL;
print("Perform the analysis with binary DNA accessibility and include weak binding sites")
source("ScriptAccessibilityRegionsWeakBinding.R")

#7. perform the analysis with non-binary DNA accessibility
profile=NULL;
bindingEnergy=NULL;
occupancy=NULL;
DNAAccessibility=NULL;
print("Perform the analysis with non-binary DNA accessibility")
source("ScriptAccessibilityProbability.R")

#8. perform the analysis with binary DNA accessibility and use the JASPAR PWM motifs
profile=NULL;
bindingEnergy=NULL;
occupancy=NULL;
DNAAccessibility=NULL;
print("Perform the analysis with binary DNA accessibility and use the JASPAR PWM motifs")
source("ScriptAccessibilityRegionsJASPAR.R")

#9. perform the genome wide analysis
profile=NULL;
bindingEnergy=NULL;
occupancy=NULL;
DNAAccessibility=NULL;
print("Perform the genome-wide analysis")
source("ScriptAccessibilityRegionsGenomeWideBestSolution.R")

