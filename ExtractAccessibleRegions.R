

if(!file.exists("objects/Dm3S5AccessibilityRegions.RData")){
	source("ProcessAccessibility.R")
}
if(file.exists("objects/Dm3S5AccessibilityRegions.RData")){
	load("objects/Dm3S5AccessibilityRegions.RData");
} else{
	stop("Could not laod DNA accessibility data");
}




#reference genome
DNASequnceSet=getSeq(Dmelanogaster,as.character=FALSE);
names(DNASequnceSet)=seqnames(Dmelanogaster);
BPFrequency=computeBPFrequency(DNASequnceSet);
#DNASequnceSet=DNASequnceSet[-union(grep("chrU", names(DNASequnceSet)),grep("Het", names(DNASequnceSet)))]# removed DNA sequences that were not assigned
DNASequnceSet=DNASequnceSet[names(DNASequnceSet)%in%c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")];# discard heterochromatin (chr2LHet, chr2RHet, chr3LHet, chr3RHet, chrXHet, chrYHet), ChrU and chrUextra (unmapped) and ChrM (mitochrondrial)


locusChromosome=c();
locusStart=c();
locusEnd=c();
fragmentLength=20000;
for(i in names(DNASequnceSet)){
	for(j in 1:ceiling(length(DNASequnceSet[[i]])/fragmentLength)){
		if(sum(dm3S5AccessibilityRegions[[i]][((j-1)*fragmentLength + 1):min(length(DNASequnceSet[[i]]),j*fragmentLength)])>0){
			locusChromosome=c(locusChromosome,i);
			locusStart = c(locusStart,((j-1)*fragmentLength + 1));
			locusEnd = c(locusEnd, min(length(DNASequnceSet[[i]]),j*fragmentLength));
		} 
	}
}
genomeWideSetAccessibleRegions=GRanges(seqnames=Rle(locusChromosome),ranges=IRanges(locusStart,locusEnd));

save(genomeWideSetAccessibleRegions,file="objects/GenomeWideSetAccessibleRegions.RData")

