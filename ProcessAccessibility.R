library("BSgenome.Dmelanogaster.UCSC.dm3")
source("functions/GenomicGeneralFunctions.R")

#reference genome
DNASequnceSet=getSeq(Dmelanogaster,as.character=FALSE);
names(DNASequnceSet)=seqnames(Dmelanogaster);
DNASequnceSet=DNASequnceSet[names(DNASequnceSet)%in%c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")];# discard heterochromatin (chr2LHet, chr2RHet, chr3LHet, chr3RHet, chrXHet, chrYHet), ChrU and chrUextra (unmapped) and ChrM (mitochrondrial)



#preprocess accessible regions
if(!file.exists("objects/Dm3S5AccessibilityRegions.RData")){
    BCDaccessibilityRaw = read.table("data/Dmel2006GeneS5Regions.tsv.gz",skip=2,header=FALSE,stringsAsFactors = FALSE);
    dm3S5AccessibilityRegions=vector("list",length(names(DNASequnceSet)));
    names(dm3S5AccessibilityRegions)=names(DNASequnceSet);
    for(i in names(DNASequnceSet)){
        dm3S5AccessibilityRegions[[i]]=rep(0,width(DNASequnceSet[i]));
    }
    chrID=match(BCDaccessibilityRaw[,2],names(DNASequnceSet));
    chrIDNotNA=which(!is.na(chrID));
    
    for(i in chrIDNotNA){
        dm3S5AccessibilityRegions[[chrID[i]]][BCDaccessibilityRaw[i,3]:BCDaccessibilityRaw[i,4]] = 1;
    }
    save(dm3S5AccessibilityRegions, file="objects/Dm3S5AccessibilityRegions.RData")
}






#preprocess probability of accessible DNA
if(!file.exists("objects/Dm3S5AccessibilityProbability.RData")){

  rep1Tab = read.table("data/Dmel2006GeneS5Rep1.bedGraph.gz",stringsAsFactors=FALSE, skip=1);
  rep2Tab = read.table("data/Dmel2006GeneS5Rep2.bedGraph.gz",stringsAsFactors=FALSE, skip=1);
  
  
  dm3S5AccessibilityReadDensity=vector("list",length(names(DNASequnceSet)));
  names(dm3S5AccessibilityReadDensity)=names(DNASequnceSet);
  for(i in names(DNASequnceSet)){
    dm3S5AccessibilityReadDensity[[i]]=rep(0,length(Dmelanogaster[[i]]));
  }
  
  chrID=match(rep1Tab[,1],names(DNASequnceSet));
  chrIDNotNA=which(!is.na(chrID));
  
  for(i in chrIDNotNA){
    chr = chrID[i];
    start = rep1Tab[i,2]-9;
    end = rep1Tab[i,3]+9;
    value=rep1Tab[i,5];
    dm3S5AccessibilityReadDensity[[chr]][start:end] = value;
  }
  
  chrID=match(rep2Tab[,1],names(DNASequnceSet));
  chrIDNotNA=which(!is.na(chrID));
  
  for(i in chrIDNotNA){
    chr = chrID[i];
    start = rep2Tab[i,2]-9;
    end = rep2Tab[i,3]+9;
    value=rep2Tab[i,5];
    dm3S5AccessibilityReadDensity[[chr]][start:end] = dm3S5AccessibilityReadDensity[[chr]][start:end]+value;
  }
  
  alpha=6.008;
  beta=0.207;
  dm3S5AccessibilityProbability=vector("list",length(dm3S5AccessibilityReadDensity));
  names(dm3S5AccessibilityProbability)=names(dm3S5AccessibilityReadDensity);
  for(index in 1:length(dm3S5AccessibilityReadDensity)){
    if(!is.null(dm3S5AccessibilityReadDensity[[index]])){
      dm3S5AccessibilityProbability[[index]]=processReadDensity((dm3S5AccessibilityReadDensity[[index]]/2), alpha, beta)
    }
  }
  save(dm3S5AccessibilityProbability, file="objects/Dm3S5AccessibilityProbability.RData");
  
}
