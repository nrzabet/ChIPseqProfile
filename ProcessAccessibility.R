library("BSgenome.Dmelanogaster.UCSC.dm3")
source("functions/GenomicGeneralFunctions.R")


#preprocess accessible regions
if(!file.exists("objects/Dm3S5AccessibilityRegions.RData")){
  BCDaccessibilityRaw = read.table("data/Dmel2006GeneS5Regions.tsv.gz", skip=2);
  dm3S5AccessibilityRegions=vector("list",length(seqnames(Dmelanogaster)));
  names(dm3S5AccessibilityRegions)=seqnames(Dmelanogaster);
  for(i in seqnames(Dmelanogaster)){
    dm3S5AccessibilityRegions[[i]]=rep(0,length(Dmelanogaster[[i]]));
  }
  for(i in 1:nrow(BCDaccessibilityRaw)){
    dm3S5AccessibilityRegions[[BCDaccessibilityRaw[i,2]]][BCDaccessibilityRaw[i,3]:BCDaccessibilityRaw[i,4]] = 1;
  }
  save(dm3S5AccessibilityRegions, file="objects/Dm3S5AccessibilityRegions.RData")
}


#preprocess probability of accessible DNA
if(!file.exists("objects/Dm3S5AccessibilityProbability.RData")){

  rep1Tab = read.table("data/Dmel2006GeneS5Rep1.bedGraph.gz",stringsAsFactors=FALSE, skip=1);
  rep2Tab = read.table("data/Dmel2006GeneS5Rep2.bedGraph.gz",stringsAsFactors=FALSE, skip=1);
  
  
  dm3S5AccessibilityReadDensity=vector("list",length(seqnames(Dmelanogaster)));
  names(dm3S5AccessibilityReadDensity)=seqnames(Dmelanogaster);
  for(i in seqnames(Dmelanogaster)){
    dm3S5AccessibilityReadDensity[[i]]=rep(0,length(Dmelanogaster[[i]]));
  }
  
  
  for(i in 2:(nrow(rep1Tab)-1)){
    chr = rep1Tab[i,1];
    start = ceiling((rep1Tab[i-1,2]+rep1Tab[i,2])/2);
    end = floor((rep1Tab[i,2]+rep1Tab[i+1,2])/2);
    value=rep1Tab[i,5];
    dm3S5AccessibilityReadDensity[[chr]][start:end] = value;
  }
  
  
  for(i in 2:(nrow(rep2Tab)-1)){
    chr = rep2Tab[i,1];
    start = ceiling((rep2Tab[i-1,2]+rep2Tab[i,2])/2);
    end = floor((rep2Tab[i,2]+rep2Tab[i+1,2])/2);
    value=rep2Tab[i,5];
    dm3S5AccessibilityReadDensity[[chr]][start:end] = dm3S5AccessibilityReadDensity[[chr]][start:end]+value;
  }
  
  alpha=6.008;
  beta=0.207;
  dm3S5AccessibilityReadDensity[["chrdmel_mitochondrion_genome"]]=NULL;
  dm3S5AccessibilityProbability=vector("list",length(dm3S5AccessibilityReadDensity));
  names(dm3S5AccessibilityProbability)=names(dm3S5AccessibilityReadDensity);
  for(index in 1:length(dm3S5AccessibilityReadDensity)){
    if(!is.null(dm3S5AccessibilityReadDensity[[index]])){
      dm3S5AccessibilityProbability[[index]]=processReadDensity((dm3S5AccessibilityReadDensity[[index]]/2), alpha, beta)
    }
  }
  save(dm3S5AccessibilityProbability, file="objects/Dm3S5AccessibilityProbability.RData")
}
