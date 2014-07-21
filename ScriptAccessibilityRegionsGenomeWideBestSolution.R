################################################################################
# LOAD FUNCTIONS
################################################################################
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm3));
suppressPackageStartupMessages(library(xtable));


source('functions/GenomicGeneralFunctions.R');
source('functions/GenomicOutputFunctions.R');
source('functions/ComputeOccupancyFunctions.R');


################################################################################
# DEFAULT PARAMETERS
################################################################################

#parameters of the model of the ChIP-seq profile; see Kaplan et al. 2011
chipMean=200; 
chipSd=200;
chipSmooth = 250;



if(!file.exists("objects/GenomeWideSetAccessibleRegions.RData")){
	source("ExtractAccessibleRegions.R");
}
if(file.exists("objects/GenomeWideSetAccessibleRegions.RData")){
	load("objects/GenomeWideSetAccessibleRegions.RData");
} else{
	stop("Could not laod the genome wide set of accessible loci");
}


names(genomeWideSetAccessibleRegions)=paste("set",1:length(genomeWideSetAccessibleRegions),sep="");



#files containing the PFM of the five TFs in raw format (a 4 row matrix, with rows: A, C, G, T and no header)
PFMFilenames = c("BCDSlx.pfm", "CADSlx.pfm", "GTSlx.pfm", "HBSlx.pfm", "KRSlx.pfm");

#TFs analysed
TF = c("BCD", "CAD","GT","HB","KR");

#threshold of the PWM
PWMThreshold=c(0.7);

#scaling factor for PWM score; see Berg and von Hippel, 1987
lambdas=c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5 ,4 ,4.5, 5);
lambdasMatrix=matrix(rep(lambdas,times=length(PFMFilenames)), ncol = length(PFMFilenames));


# number of molecules bound to the DNA
boundMolecules = c(1, 10, 20, 50, 100, 200, 500,1000,2000, 5000,10000,20000,50000, 100000, 200000, 500000, 1000000);
boundMoleculesMatrix = matrix(rep(boundMolecules,times=length(PFMFilenames)), ncol = length(PFMFilenames));


#load already computed occupancy data for the 21 loci
load(file=paste("objects/occupancyAllTFsAllPositiveAccessibilityRegions.RData",sep=""));
lambdasMatrixNew = matrix(rep(0,times=length(PFMFilenames)), ncol = length(PFMFilenames));
boundMoleculesMatrixNew = matrix(rep(0,times=length(PFMFilenames)), ncol = length(PFMFilenames));
for(i in 1:length(PWMThreshold)){
	for(TFid in 1: length(TF)){
		optimalParameters = getOptimalSetOfParameters(occupancy[[i]], c(TF[TFid]), lambdasMatrix[,TFid], boundMoleculesMatrix[,TFid]);
		l=optimalParameters[1,1];
		m=optimalParameters[1,2];
		lambdasMatrixNew[1,TFid]=lambdasMatrix[l,TFid];
		boundMoleculesMatrixNew[1,TFid]=boundMoleculesMatrix[m,TFid];
	}
}

lambdasMatrix=lambdasMatrixNew;
boundMoleculesMatrix=boundMoleculesMatrixNew;

# the size of the TFs. If NULL, this is approximated by the length of the PWM
TFSize = NULL;

#files containing the ChIP-seq profiles
profileFiles= list("BCD"="data/BCD11.wig.gz","CAD"="data/CAD11.wig.gz","GT"="data/GT21.wig.gz","HB"="data/HB12.wig.gz","KR"="data/KR11.wig.gz");



#reference genome
DNASequnceSet=getSeq(Dmelanogaster,as.character=FALSE);
names(DNASequnceSet)=seqnames(Dmelanogaster);

#

################################################################################
# LOAD GENE REFERENCE
################################################################################
loadGeneRef=FALSE;
if(sum(ls(envir = .GlobalEnv) == "geneRef") == 0){
	loadGeneRef=TRUE;
} else if(is.null(geneRef)){
	loadGeneRef=TRUE;
}

if(loadGeneRef){
	print("(Re)loading D.melanogaster gene reference");
	if(file.exists("objects/Dmel3GeneRef.RData")){
            load("objects/Dmel3GeneRef.RData");
        } else{
		geneRefFilename = "data/DmelGenesR526.tsv.gz";
		geneRef = extractGeneList(geneRefFilename);
		save(geneRef, file="objects/Dmel3GeneRef.RData");
        }

} else{
	print("Using previously loaded D.melanogaster gene reference")
}





################################################################################
# LOAD CHIP-SEQ DATA
################################################################################
print("ChIP-seq data")
loadProfiles=FALSE;
if(sum(ls(envir = .GlobalEnv) == "profile") == 0){
	loadProfiles=TRUE;
} else if(is.null(profile)){
	loadProfiles=TRUE;
}
if(loadProfiles){
	if(file.exists("objects/ChIPSeqProfilesGenomeWide.RData") & file.exists("objects/ChIPSeqProfilesBackground.RData") & file.exists("objects/ChIPSeqProfilesMaxSignal.RData")){
        load("objects/ChIPSeqProfilesGenomeWide.RData");
        load("objects/ChIPSeqProfilesBackground.RData");
        load("objects/ChIPSeqProfilesMaxSignal.RData");
    } else{
		profile=vector("list",length(profileFiles));
		names(profile)=names(profileFiles);
		backgroundSignal=vector("list",length(profileFiles));
		names(backgroundSignal)=names(profileFiles);
		maxSignal=vector("list",length(profileFiles));
		names(maxSignal)=names(profileFiles);
		for(i in 1:length(profileFiles)){
			profile[[i]] = read.table(profileFiles[[i]], skip=2);
			profile[[i]][,4] = as.numeric(profile[[i]][,4]);
			backgroundSignal[[i]] = mean(profile[[i]][,4]);
			maxSignal[[i]] = max(profile[[i]][,4]);
			profile[[i]] =  extractOccupancyDataAtLoci(profile=profile[[i]], setSequence=genomeWideSetAccessibleRegions, maxSignal=maxSignal[[i]], removeBackground=0, chipSmooth=chipSmooth);
		}
        save(profile, file="objects/ChIPSeqProfilesGenomeWide.RData");
		save(backgroundSignal, file="objects/ChIPSeqProfilesBackground.RData");
		save(maxSignal, file="objects/ChIPSeqProfilesMaxSignal.RData");
	}
} else{
	print("Using previously loaded ChIP-seq profiles")
}



################################################################################
# LOAD DNA ACCESSIBILITY
################################################################################
loadDNAAccessibility=FALSE;
if(sum(ls(envir = .GlobalEnv) == "dm3S5AccessibilityRegions") == 0){
	loadDNAAccessibility=TRUE;
} else if(is.null(dm3S5AccessibilityRegions)){
	loadDNAAccessibility=TRUE;
}

if(loadDNAAccessibility){
	print("(Re)loading DNA accessibility data");
	if(!file.exists("objects/Dm3S5AccessibilityRegions.RData")){
		source("ProcessAccessibility.R")
	}
	if(file.exists("objects/Dm3S5AccessibilityRegions.RData")){
		load("objects/Dm3S5AccessibilityRegions.RData");
	} else{
		stop("Could not laod DNA accessibility data");
	}
} else{
	print("Using previously loaded DNA accessibility data")
}

################################################################################
# COMPUTE PWM SCORES
################################################################################
toComputePWMScores=FALSE;
if(sum(ls(envir = .GlobalEnv) == "bindingEnergy") == 0){
	toComputePWMScores=TRUE;
} else if(is.null(bindingEnergy)){
	toComputePWMScores=TRUE;
}

if(toComputePWMScores){
	print("(Re)computing PWM scores")
	
	if(file.exists("objects/BindingEnergyGenomeWideAccessibilityRegions.RData")){
		load("objects/BindingEnergyGenomeWideAccessibilityRegions.RData");
        } else{
		bindingEnergy=vector("list", length(PWMThreshold));
		for(i in 1:length(PWMThreshold)){
			bindingEnergy[[i]]=vector("list", nrow(lambdasMatrix));
			for(j in 1:nrow(lambdasMatrix)){
				bindingEnergy[[i]][[j]] = computePWMScores(directory="data", PFMFilename=PFMFilenames, PFMFormat="raw", DNASequenceSet=DNASequnceSet, setSequence=genomeWideSetAccessibleRegions,
									   TF=TF, DNAAccessibility=dm3S5AccessibilityRegions, PWMPseudocount=1,
									   PWMUseNaturalLog=FALSE, PWMNoOfSites=NULL, PWMThreshold=rep(PWMThreshold[i],times=length(TF)),lambda=lambdasMatrix[j,],
									   strand="+-", strandRule="max");
			}
		}
		save(bindingEnergy, file=paste("objects/BindingEnergyGenomeWideAccessibilityRegions.RData",sep=""));

	}	
} else{
	print("Using previously computed PWM scores");
}



################################################################################
# EXTRACT DNA ACCESSIBILITY AT LOCI
################################################################################
extractDNAAccessibilityLoci=FALSE;
if(sum(ls(envir = .GlobalEnv) == "DNAAccessibility") == 0){
	extractDNAAccessibilityLoci=TRUE;
} else if(is.null(DNAAccessibility)){
	extractDNAAccessibilityLoci=TRUE;
} else if(length(DNAAccessibility)!=length(genomeWideSetAccessibleRegions)){
	extractDNAAccessibilityLoci=TRUE;
}

if(extractDNAAccessibilityLoci){
	print("Extracting DNA accessibility at loci")
	DNAAccessibility=extractDNAAccesibilityAtLoci(dm3S5AccessibilityRegions,genomeWideSetAccessibleRegions);
} else{
	print("Using previously extracted DNA accessibility at loci");
}



################################################################################
# COMPUTE OCCUPANCY
################################################################################
toComputeOccupancy=FALSE;
if(sum(ls(envir = .GlobalEnv) == "occupancyGenomeWide") == 0){
	toComputeOccupancy=TRUE;
} else if(is.null(occupancyGenomeWide)){
	toComputeOccupancy=TRUE;
}

if(toComputeOccupancy){
	print("(Re)computing genomic occupancy of TFs")
	if(file.exists("objects/occupancyAllTFsAllPositiveAccessibilityRegionsGenomeWide.RData")){
		load("objects/occupancyAllTFsAllPositiveAccessibilityRegionsGenomeWide.RData");
        } else{
		occupancyGenomeWide=vector("list", length(PWMThreshold));
		for(i in 1:length(PWMThreshold)){
			occupancyGenomeWide[[i]]=vector("list", nrow(lambdasMatrix));
			for(j in 1:nrow(lambdasMatrix)){
				occupancyGenomeWide[[i]][[j]]=vector("list", nrow(boundMoleculesMatrix));
				for(k in 1:nrow(boundMoleculesMatrix)){
				print(paste("PWMthreshold = ",(formatC(round(100*PWMThreshold[i]), width = 3, format = "d", flag = "0")),"; lambda = ",(formatC(round(100*lambdasMatrix[j,]), width = 3, format = "d", flag = "0")),"; abundance=", boundMoleculesMatrix[k,], sep=""));
				occupancyGenomeWide[[i]][[j]][[k]] = computeOccupancy(TF=TF, PWMScore=bindingEnergy[[i]][[j]]$PWMScore, indexPWMThresholded=bindingEnergy[[i]][[j]]$indexPWMThresholded,
									    DNALength=bindingEnergy[[i]][[j]]$DNALength, averageExpPWMScore=bindingEnergy[[i]][[j]]$averageExpPWMScore, DNAAccessibility=DNAAccessibility,
									    setSequence=genomeWideSetAccessibleRegions, lambda=lambdasMatrix[j,], boundMolecules=boundMoleculesMatrix[k,], norm=TRUE,
									    outputBoundMoleculesOccupancy=TRUE, backgroundSignal=unlist(backgroundSignal), maxSignal=unlist(maxSignal), outputChIPseqProfile=TRUE,
									    chipMean= chipMean, chipSd=chipSd, chipSmooth = chipSmooth, removeBackground=0, profile=profile);
				}	
			}
		}
		save(occupancyGenomeWide, file=paste("objects/occupancyAllTFsAllPositiveAccessibilityRegionsGenomeWide.RData",sep=""));
		
	} 
} else{
	print("Using previously computed genomic occupancy of TFs")
}


################################################################################
# PLOT STATISTICS
################################################################################
for(imgType in c("pdf")){
	GenomeWideNoOfRegions=plotStatistics(occupancyGenomeWide[[1]][[1]][[1]], TF, directory="img/statistics", plotFilename="BestSolutionAllAccessibleRegions", imageType=imgType, threshold=c(1.0,0.5),backgroundSignal,profile);
}




################################################################################
# NUMBER OF DNA SEGMENTS WITH A CHIP-SEQ SIGNAL HIGHER THAN THE THRESHOLD
################################################################################
tableCaption="\\emph{The number of DNA segments with a ChIP-seq signal higher than the threshold}. We considered two thresholds $K=1.0$ and $K=0.5$. The mean ChIP-seq signal of the $20\\ Kbp$ segment needs to be higher than $K\\times B$.";	
tableLatex=xtable(GenomeWideNoOfRegions, label ="tab:GenomeWideNoOfRegions", caption =tableCaption,align="|l|r|r|r|r|r|", digit=c(0));
print.xtable(tableLatex, sanitize.text.function = function(x) x,file="tables/GenomeWideNoOfRegions.tex",hline.after=-1:nrow(GenomeWideNoOfRegions),floating.environment="table");


################################################################################
# ChIP-seq PROFILE STATISTICS
################################################################################
ChIPseqProfileStatistics = matrix(0,ncol=length(TF),nrow=2);
colnames(ChIPseqProfileStatistics)=TF;
row.names(ChIPseqProfileStatistics)=c("$B$", "$M$");
ChIPseqProfileStatistics[1,]=unlist(backgroundSignal);
ChIPseqProfileStatistics[2,]=unlist(maxSignal);
tableCaption="\\emph{Mean and maximum of the ChIP-seq signal for the five TFs}.";	
tableLatex=xtable(ChIPseqProfileStatistics, label ="tab:ChIPseqProfileStatistics", caption =tableCaption,align="|l|r|r|r|r|r|");
print.xtable(tableLatex, sanitize.text.function = function(x) x,file="tables/ChIPseqProfileStatistics.tex",hline.after=-1:nrow(ChIPseqProfileStatistics),floating.environment="table");


