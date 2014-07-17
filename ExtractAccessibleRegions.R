

if(!file.exists("objects/Dm3S5AccessibilityRegions.RData")){
	source("ProcessAccessibility.R")
}
if(file.exists("objects/Dm3S5AccessibilityRegions.RData")){
	load("objects/Dm3S5AccessibilityRegions.RData");
} else{
	stop("Could not laod DNA accessibility data");
}






locusChromosome=c();
locusStart=c();
locusEnd=c();
fragmentLength=20000;
for(i in seqnames(Dmelanogaster)){
	for(j in 1:ceiling(length(Dmelanogaster[[i]])/fragmentLength)){
		if(sum(dm3S5AccessibilityRegions[[i]][((j-1)*fragmentLength + 1):min(length(Dmelanogaster[[i]]),j*fragmentLength)])>0){
			locusChromosome=c(locusChromosome,i);
			locusStart = c(locusStart,((j-1)*fragmentLength + 1));
			locusEnd = c(locusEnd, min(length(Dmelanogaster[[i]]),j*fragmentLength));
		} 
	}
}
genomeWideSetAccessibleRegions=GRanges(seqnames=Rle(locusChromosome),ranges=IRanges(locusStart,locusEnd));

save(genomeWideSetAccessibleRegions,file="objects/GenomeWideSetAccessibleRegions.RData")

