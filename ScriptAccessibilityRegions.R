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

#the 21 loci of interest
allSetPositive=GRanges(seqnames=Rle(c("chr3L","chr3R","chr2L","chr3L","chr3R","chr3L","chr2L","chr3L","chr2R","chr2L","chrX","chr3R","chr3R","chrX","chr3R","chr2R","chr2L","chrX","chr3R","chr3R","chrX")), 
ranges=IRanges(c(21461001,19011001,3820001,20683260,169001,14165001,12077501,8656154,5860693,20767501,8518001,670001,2682501,2317878,4513501,21103924,3603001,20548001,24396001,26672001,18193001),
c(21477000,19024000,3840000,20695259,181000,14179000,12095500 ,8682153 ,5876692 ,20786500 ,8550000 ,696000 ,2696500 ,2330877 ,4531500 ,21118923 ,3613000 ,20570000 ,24420000 ,26684000 ,18208000)));
names(allSetPositive)=c("croc","cnc","slp","kni","hkb","D","prd","H","eve","cad","oc","opa", "ftz", "gt", "hb", "Kr", "odd", "run", "fkh", "tll", "os");


eveLocus=GRanges(seqnames=Rle(c("chr2R")), ranges=IRanges(c(5860693),c(5876692)));
names(eveLocus)=c("eve");

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
	if(file.exists("objects/ChIPSeqProfilesAtLoci.RData") & file.exists("objects/ChIPSeqProfilesBackground.RData") & file.exists("objects/ChIPSeqProfilesMaxSignal.RData")){
            load("objects/ChIPSeqProfilesAtLoci.RData");
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
		  profile[[i]] =  extractOccupancyDataAtLoci(profile=profile[[i]], setSequence=allSetPositive, maxSignal=maxSignal[[i]], removeBackground=0, chipSmooth=chipSmooth);
		}
	  save(profile, file="objects/ChIPSeqProfilesAtLoci.RData");
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
	
	if(file.exists("objects/BindingEnergyAtLociAccessibilityRegions.RData")){
		load("objects/BindingEnergyAtLociAccessibilityRegions.RData");
        } else{
		bindingEnergy=vector("list", length(PWMThreshold));
		for(i in 1:length(PWMThreshold)){
			bindingEnergy[[i]]=vector("list", nrow(lambdasMatrix));
			for(j in 1:nrow(lambdasMatrix)){
				bindingEnergy[[i]][[j]] = computePWMScores(directory="data", PFMFilename=PFMFilenames, PFMFormat="raw", DNASequenceSet=DNASequnceSet, setSequence=allSetPositive,
									   TF=TF, DNAAccessibility=dm3S5AccessibilityRegions, PWMPseudocount=1,
									   PWMUseNaturalLog=FALSE, PWMNoOfSites=NULL, PWMThreshold=rep(PWMThreshold[i],times=length(TF)),lambda=lambdasMatrix[j,],
									   strand="+-", strandRule="max");
			}
		}
		save(bindingEnergy, file=paste("objects/BindingEnergyAtLociAccessibilityRegions.RData",sep=""));

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
} else if(length(DNAAccessibility)!=length(allSetPositive)){
	extractDNAAccessibilityLoci=TRUE;
}

if(extractDNAAccessibilityLoci){
	print("Extracting DNA accessibility at loci")
	DNAAccessibility=extractDNAAccesibilityAtLoci(dm3S5AccessibilityRegions,allSetPositive);
} else{
	print("Using previously extracted DNA accessibility at loci");
}



################################################################################
# COMPUTE OCCUPANCY
################################################################################
toComputeOccupancy=FALSE;
if(sum(ls(envir = .GlobalEnv) == "occupancy") == 0){
	toComputeOccupancy=TRUE;
} else if(is.null(occupancy)){
	toComputeOccupancy=TRUE;
}

if(toComputeOccupancy){
	print("(Re)computing genomic occupancy of TFs")
	if(file.exists("objects/occupancyAllTFsAllPositiveAccessibilityRegions.RData")){
		load("objects/occupancyAllTFsAllPositiveAccessibilityRegions.RData");
        } else{
		occupancy=vector("list", length(PWMThreshold));
		for(i in 1:length(PWMThreshold)){
			occupancy[[i]]=vector("list", nrow(lambdasMatrix));
			for(j in 1:nrow(lambdasMatrix)){
				occupancy[[i]][[j]]=vector("list", nrow(boundMoleculesMatrix));
				for(k in 1:nrow(boundMoleculesMatrix)){
					print(paste("PWMthreshold = ",(formatC(round(100*PWMThreshold[i]), width = 3, format = "d", flag = "0")),"; lambda = ",(formatC(round(100*lambdasMatrix[j,]), width = 3, format = "d", flag = "0")),"; abundance=", boundMoleculesMatrix[k,], sep=""));
					occupancy[[i]][[j]][[k]] = computeOccupancy(TF=TF, PWMScore=bindingEnergy[[i]][[j]]$PWMScore, indexPWMThresholded=bindingEnergy[[i]][[j]]$indexPWMThresholded,
							    DNALength=bindingEnergy[[i]][[j]]$DNALength, averageExpPWMScore=bindingEnergy[[i]][[j]]$averageExpPWMScore, DNAAccessibility=DNAAccessibility,
							    setSequence=allSetPositive, lambda=lambdasMatrix[j,], boundMolecules=boundMoleculesMatrix[k,], norm=TRUE,
							    outputBoundMoleculesOccupancy=TRUE, backgroundSignal=unlist(backgroundSignal), maxSignal=unlist(maxSignal), outputChIPseqProfile=TRUE,
							    chipMean= chipMean, chipSd=chipSd, chipSmooth = chipSmooth, removeBackground=0, profile=profile);
				}	
			}
		}
		save(occupancy, file=paste("objects/occupancyAllTFsAllPositiveAccessibilityRegions.RData",sep=""));
		
	} 
} else{
	print("Using previously computed genomic occupancy of TFs")
}


################################################################################
# PLOT HEATMAPS
################################################################################
for(imageType in c("pdf")){
	for(i in 1:length(PWMThreshold)){
		txtCase = paste("PWMthreshold",(formatC(round(100*PWMThreshold[i]), width = 3, format = "d", flag = "0")),sep="");	
		plotOccupancyModelQualityHeatmaps(occupancy=occupancy[[i]], TF=c("BCD","CAD"), lambdas=lambdasMatrix[,1], boundMolecules=boundMoleculesMatrix[,1], directory="img/heatmap",
						  plotFilename=paste("BCDCADAccessibiltiyRegions",txtCase,sep=""), imageType=imageType, contour=TRUE);
		plotOccupancyModelQualityHeatmaps(occupancy=occupancy[[i]], TF=c("GT","HB","KR"), lambdas=lambdasMatrix[,1], boundMolecules=boundMoleculesMatrix[,1], directory="img/heatmap",
						  plotFilename=paste("GTHBKRAccessibiltiyRegions",txtCase,sep=""), imageType=imageType, contour=TRUE);
	}
}


################################################################################
# OPTIMAL SET OF PARAMETERS
################################################################################
optimalSetOfParameters = matrix(0,nrow=length(TF),ncol=5);
row.names(optimalSetOfParameters)=TF;
colnames(optimalSetOfParameters)=c("N", "$\\lambda$", "$\\theta$","$MSE$", "$\\rho$");
for(i in 1:length(PWMThreshold)){
	for(TFid in 1: length(TF)){
		optimalParameters = getOptimalSetOfParameters(occupancy[[i]], c(TF[TFid]), lambdasMatrix[,TFid], boundMoleculesMatrix[,TFid]);
		l=optimalParameters[1,1];
		m=optimalParameters[1,2];
		
		optimalSetOfParameters[TFid,1]=boundMoleculesMatrix[m,TFid];
		optimalSetOfParameters[TFid,2]=format(round(lambdasMatrix[l,TFid], 2), nsmall = 2);
		optimalSetOfParameters[TFid,3]=format(round(occupancy[[i]][[l]][[m]]$meanTheta[[TF[TFid]]], 2), nsmall = 2);
		
		optimalMSE=format(round(occupancy[[i]][[l]][[m]]$meanMSE[[TF[TFid]]], 2), nsmall = 2);
		optimalCorrelation=format(round(occupancy[[i]][[l]][[m]]$meanCorrelation[[TF[TFid]]], 2), nsmall = 2);
		
		optimalParameters = getOptimalSetOfParameters(occupancy[[i]], c(TF[TFid]), lambdasMatrix[,TFid], boundMoleculesMatrix[,TFid],parameter="MSE");
		l=optimalParameters[1,1];
		m=optimalParameters[1,2];
		minMSE=format(round(occupancy[[i]][[l]][[m]]$meanMSE[[TF[TFid]]], 2), nsmall = 2);

		optimalParameters = getOptimalSetOfParameters(occupancy[[i]], c(TF[TFid]), lambdasMatrix[,TFid], boundMoleculesMatrix[,TFid],parameter="correlation");
		l=optimalParameters[1,1];
		m=optimalParameters[1,2];
		maxCorrelation=format(round(occupancy[[i]][[l]][[m]]$meanCorrelation[[TF[TFid]]], 2), nsmall = 2);

		optimalSetOfParameters[TFid,4]=paste(optimalMSE," (",minMSE,")",sep="");
		optimalSetOfParameters[TFid,5]=paste(optimalCorrelation," (",maxCorrelation,")",sep="");

	}
}
	
tableCaption="\\emph{Set of parameters that minimises the difference between the ChIP-seq profile and the analytical model which includes DNA accessibility}. The accessibility of any site can be either $0$ or $1$ depending on whether the site is accessible or not. We also listed the values for the mean square error ($MSE$) and correlation ($\\rho$). The values in the parentheses represent the minimum mean square error and the maximum correlation. We considered only the sites that have a PWM score higher than $70\\%$ of the distance between the lowest and the highest score.";	
tableLatex=xtable(optimalSetOfParameters, label ="tab:paramsModelAccRegions", caption =tableCaption,align="|l|r|r|r|r|r|");
print.xtable(tableLatex, sanitize.text.function = function(x) x,file="tables/OptimalSetOfParametersTableAccessibilityRegions.tex",hline.after=-1:nrow(optimalSetOfParameters),floating.environment="table")




################################################################################
# PLOT PROFILES
################################################################################
for(imageType in c("pdf")){
	for(i in 1:length(PWMThreshold)){
		for(TFid in 1: length(TF)){
			optimalParameters = getOptimalSetOfParameters(occupancy[[i]], c(TF[TFid]), lambdasMatrix[,TFid], boundMoleculesMatrix[,TFid]);
			l=optimalParameters[1,1];
			m=optimalParameters[1,2];
			txtCase = paste("PWMThreshold",(formatC(round(100*PWMThreshold[i]), width = 3, format = "d", flag = "0")),
					"lambda",(formatC(round(100*lambdasMatrix[l,TFid]), width = 3, format = "d", flag = "0")), sep="");
			plotOccupancyProfile(occupancy=occupancy[[i]][[l]][[m]]$occupancyAbundanceChIP[[TFid]], PWMScore=bindingEnergy[[i]][[l]]$PWMScore[[TFid]], lambda=lambdasMatrix[l,TFid],
					     maxPWMScore=bindingEnergy[[i]][[l]]$maxPWMScore[[TFid]], DNAAccessibility=DNAAccessibility, directory="img/profile",
					     plotFilename=paste(TF[TFid],"AccessibilityRegions",txtCase,sep=""), profile=profile[[TFid]], imageType=imageType, setSequence=allSetPositive,
					     TF=TF[TFid], boundMolecules=boundMoleculesMatrix[m,TFid], outputPWMScore=TRUE, outputOccupancyAbundance=TRUE, geneRef=geneRef, stepSize=20);
		}		
	}
}



################################################################################
# PLOT PROFILES AT EVE STRIPE ENHANCER
################################################################################
geneRefNew=geneRef;
enhancerFile = "data/EveStripeEnhancer.bed.gz";
bufferEnhancer = extractEnhancers(enhancerFile);
geneRefNew$enhancer=bufferEnhancer;
for(imageType in c("pdf")){
	for(i in 1:length(PWMThreshold)){
		for(TFid in which(TF=="BCD")){
			optimalParameters = getOptimalSetOfParameters(occupancy[[i]], c(TF[TFid]), lambdasMatrix[,TFid], boundMoleculesMatrix[,TFid]);
			l=optimalParameters[1,1];
			for(m in which(boundMoleculesMatrix[,TFid]%in%c(20000,1000,50))){
				txtCase = paste("PWMThreshold",(formatC(round(100*PWMThreshold[i]), width = 3, format = "d", flag = "0")),
						"lambda",(formatC(round(100*lambdasMatrix[l,TFid]), width = 3, format = "d", flag = "0")), sep="");
				plotOccupancyProfile(occupancy=occupancy[[i]][[l]][[m]]$occupancyAbundanceChIP[[TFid]], PWMScore=bindingEnergy[[i]][[l]]$PWMScore[[TFid]], lambda=lambdasMatrix[l,TFid],
							maxPWMScore=bindingEnergy[[i]][[l]]$maxPWMScore[[TFid]], DNAAccessibility=DNAAccessibility, directory="img/profile",
							plotFilename=paste(TF[TFid],"AccessibilityRegionsEveEnhancer",txtCase,sep=""), profile=NULL, imageType=imageType, setSequence=eveLocus,
							TF=TF[TFid], boundMolecules=boundMoleculesMatrix[m,TFid], outputPWMScore=TRUE, outputOccupancyAbundance=TRUE, geneRef=geneRefNew, stepSize=20);
			}
		}		
	}
}
