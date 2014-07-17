

if(!require(BSgenome.Dmelanogaster.UCSC.dm3)){
    source("http://bioconductor.org/biocLite.R");
    biocLite("BSgenome.Dmelanogaster.UCSC.dm3");
}

if(!require(xtable)){
    install.packages("xtable");
}
