#library("Biostrings");
#library("GenomicRanges");
#

### computes the PWM matrices for a list of TFs and then scores the entire genome and computes the average waiting time
### directory: the directory containing the PFM files
### PFMFilename: an array containing the files taht store the PFM for each TFs. It uses the following formats: RAW, TRANSFAC, JASPAR or constructs the PFM from a set of sequences. For the correct formats see http://www.benoslab.pitt.edu/stamp/help.html#input
### PFMFormat="raw": the format of the PWM: raw, transfac, JASPAR, sequences (list of binding sequences)
### DNASequenceSet: a DNAStringSet containing the DNA sequences
### setSequence: a list of GRanges where the occupancy will be computed. If NULL, the occupancy is computed genome wide
### TF=NULL: the names of the TFs
### DNAAccessibility=NULL: the file containing the DNA accessibility data. If NULL, all DNA is assumed to be accessible.
### PWMPseudocount=1: the pseudo count value used when computing the PWM
### PWMUseNaturalLog=FALSE: this is a logical value indicating whether to use natural logharithm when computing the PWM. If FALSE, the PWM is computed using logharithm in base 2.
### PWMNoOfSites=NULL: the number of sites used to construc the PWM. If NULL, it uses the PFM to compute this.
### PWMThreshold=NULL: the threshold for the PWM. If NULL, it considers all sites.
### lambda=1: the scalling factor between the PWM and the binding energy E=w*(1/lambda).
### strand="+": the strand used to compute the PWM score; can be "+", "-", "+-"
### strandRule="max": the rule used to combine the PWM score on two strands; can be "max", "sum" or "mean"
### return: a list consisting of:
###         DNALength: the length of the genome
###         PWMScore: a list (for each TF) consisting of a list (for each locus in the setSequence list) with the PWM score
###         maxPWMScore: a list consisting of the maximum PWM score for each TF
###         minPWMScore: a list consisting of the minimum PWM score for each TF
###         indexPWMThresholded: a list (for each TF) consisting of a list (for each locus in the setSequence list) with the positions where the PWM score is higher than the PWMthreshold
###         averageExpPWMScore: a list consisting of the genome-wide computed average of the DNA accessibility data multiplied with the exponential of the binding energy (PWM score multiplied to the 1/lambda) for each TF
###         PWM: a list consisting of the PFM matrix for each TF
###         PFM: a list consisting of the PWM matrix for each TF
computePWMScores <- function(directory, PFMFilename, PFMFormat="raw", DNASequenceSet, setSequence,TF=NULL, DNAAccessibility=NULL, PWMPseudocount=1, PWMUseNaturalLog=FALSE, PWMNoOfSites=NULL, PWMThreshold=NULL, lambda=1, strand="+", strandRule="max"){
	
	#the path to the PFM files
	path="";
	if(!is.null(directory) & nchar(directory) > 0){
		path=paste(directory,"/", sep="");
	}
	
	#extract the PFM from a file
	print("extract the PFM");
	PFM = vector("list",length(PFMFilename));
	if(!is.null(TF) & length(TF)==length(PFMFilename)){
		names(PFM)=TF;
	}
	
	for(i in 1:length(PFMFilename)){
		PFM[[i]] = parsePFM(paste(path,PFMFilename[i],sep=""), format=PFMFormat);
		if(is.null(PFM[[i]])){
			cat("Could not parse the PFM from file: ",PFMFilename, file=paste(path,"status.err",sep="",append=TRUE),sep="\n");
			stop(paste("Could not parse the PFM from file: ",PFMFilename));

		}
	}

	
	if(!is.null(DNASequenceSet)){
		#base pair frequency in the genome
		BPFrequency=computeBPFrequency(DNASequenceSet);
		
		#construct the PWMs
		PWM = vector("list",length(PFMFilename));
		if(!is.null(TF) & length(TF)==length(PFMFilename)){
			names(PWM)=TF;
		}
		for(i in 1:length(PFMFilename)){
			PWM[[i]] = computePWM(PFM=PFM[[i]],pseudocount=PWMPseudocount,BPFrequency= BPFrequency, naturalLog=PWMUseNaturalLog,noOfSites=PWMNoOfSites);
		}
		
		print("score the genome with the PWM");
		#score the genome with the new PWM
		DNASequenceScoreSet = vector("list",length(PFMFilename));
		if(!is.null(TF) & length(TF)==length(PFMFilename)){
			names(DNASequenceScoreSet)=TF;
		}
		for(i in 1:length(PFMFilename)){
			DNASequenceScoreSet[[i]]=scoreDNAStringSet(PWM[[i]],DNASequenceSet,strand, strandRule);
		}
		
		
		
		# compute the number of possible site, which is given by the length of the DNA 
		DNALength=unlist(lapply(DNASequenceScoreSet[[1]],length));
		DNALengthTotal = sum(DNALength);
		

		#compute mean waiting time, max PWM score and min PWM score
		print("compute the mean waiting time");
		averageExpPWMScore = vector("list",length(PFMFilename));
		maxPWMScore = vector("list",length(PFMFilename));
		minPWMScore = vector("list",length(PFMFilename));
		if(!is.null(TF) & length(TF)==length(PFMFilename)){
			names(averageExpPWMScore)=TF;
			names(maxPWMScore)=TF;
			names(minPWMScore)=TF;
		}
		
		for(i in 1:length(PFMFilename)){
			if(!is.null(DNAAccessibility) & is.list(DNAAccessibility) & length(DNASequenceScoreSet[[i]])==length(DNAAccessibility)){
				averageExpPWMScoreLocal=c();
				for(j in 1:length(DNASequenceScoreSet[[i]])){
					minLength = min(length(DNASequenceScoreSet[[i]][[j]]),length(DNAAccessibility[[j]]));	
					indexes=1:minLength;
					averageExpPWMScoreLocal = c(averageExpPWMScoreLocal,mean(DNAAccessibility[[j]][indexes]*exp((1/lambda[i])*DNASequenceScoreSet[[i]][[j]][indexes])));	
				}
			} else{
				averageExpPWMScoreLocal = unlist(lapply(lapply(lapply(DNASequenceScoreSet[[i]],"*", (1/lambda[i])),exp),mean));
			}
			
			
			averageExpPWMScore[[i]] = mean(averageExpPWMScoreLocal*(DNALength/DNALengthTotal));
			maxPWMScore[[i]] = max(unlist(lapply(DNASequenceScoreSet[[i]],max)));
			minPWMScore[[i]] = min(unlist(lapply(DNASequenceScoreSet[[i]],min)));
			
		}
		

			
		if(!is.null(setSequence) & length(setSequence)>0){
			print("extract the PWM score at the loci");
			#extract the PWM scores only at the needed positions
			DNASequenceScoreSetLocal = vector("list",length(PFMFilename));
			if(!is.null(TF) & length(TF)==length(PFMFilename)){
				names(DNASequenceScoreSetLocal)=TF;
			}
			for(i in 1:length(PFMFilename)){
				DNASequenceScoreSetBuffer = vector("list",length(setSequence));				
				names(DNASequenceScoreSetBuffer)=names(setSequence);
				
				for(gr in 1:length(setSequence)){
					DNASequenceScoreSetBuffer[[gr]] = DNASequenceScoreSet[[i]][[as.vector(seqnames(setSequence))[gr]]][as.vector(ranges(setSequence)@start)[gr]:as.vector(end(ranges(setSequence)))[gr]];	
				}
				DNASequenceScoreSetBuffer[[gr]][which(is.na(DNASequenceScoreSetBuffer[[gr]]))]=minPWMScore[[i]];
				DNASequenceScoreSetLocal[[i]]=DNASequenceScoreSetBuffer;
			}
		} else{
			DNASequenceScoreSetLocal=DNASequenceScoreSet;
			DNASequenceScoreSet = NULL;
		}
		
		#get the positions of the sites with a PWM score higher than the threshold
		indexPWMThresholded = vector("list",length(PFMFilename));
		names(indexPWMThresholded)=TF;
		for(i in 1:length(PFMFilename)){
			PWMThresholdLocal=minPWMScore[[i]] + PWMThreshold[i]*(maxPWMScore[[i]]-minPWMScore[[i]]);
			indexPWMThresholded[[i]]  = getIndexOfPWMThresholded(DNASequenceScoreSetLocal[[i]],PWMThresholdLocal);
		}

		result=list("DNALength"=DNALengthTotal,"PWMScore"=DNASequenceScoreSetLocal, "maxPWMScore"=maxPWMScore, "minPWMScore"=minPWMScore, "indexPWMThresholded"=indexPWMThresholded, "averageExpPWMScore"=averageExpPWMScore, "PWM"=PWM, "PFM"=PFM);
	}
	
	
	return(result);
}

### generates a list of all PWM scores that are above a certain threshold 
### PWMScore: a list a PWM scores vectors for each DNA segments
### PWMThreshold: the PWM score threshold
### minPWMScore: the min PWM score
### maxPWMScore: the max PWM score
### return: a list consisting of a vector with the indexes of the PWM scores higher than a threshold for each chromosome (or locus)
getIndexOfPWMThresholded <- function(PWMScore, PWMThreshold){
	indexes = vector("list",length(PWMScore));
	names(indexes)=names(PWMScore);
	for(i in 1:length(PWMScore)){
		indexes[[i]] = which(PWMScore[[i]]>PWMThreshold);
	} 
	return(indexes);
}

### extracts the DNA accessbility data for the specified loci
### DNAAccessibility: the genome-wide DNA accessibility data
### setSequence: a GenomicRanges object consisting of all the loci of interest
### setSequenceName: the names of all the loci of interest
### return: a list containing the accessibility data at each loci in the form of a vector. 
extractDNAAccesibilityAtLoci <-function(DNAAccessibility,setSequence){
	DNAAccessibilityBuffer=NULL;
	if(!is.null(DNAAccessibility)){
		#extract the DNA accessibility only at the loci
		if(!is.null(setSequence) & length(setSequence)>0){
			DNAAccessibilityBuffer = vector("list",length(setSequence));
			names(DNAAccessibilityBuffer)=names(setSequence);
		
			for(gr in 1:length(setSequence)){
				DNAAccessibilityBuffer[[gr]] = DNAAccessibility[[as.vector(seqnames(setSequence))[gr]]][as.vector(ranges(setSequence)@start)[gr]:as.vector(end(ranges(setSequence)))[gr]];	
			}	
		}
	}
	return(DNAAccessibilityBuffer);
}

### Compute the occupancy given a DNA sequence and a TF motif
### directory: the directory containing the data files
### TF: the names of the TFs
### PWMScore: a list of vectors with the PWM score at each position
### indexPWMThresholded: a list of vectors with the site with a PWM score higher than a threshold
### DNALength: the length of the genome
### averageExpPWMScore: the average of the exponential of the PWM score for each TF
### DNAAccessibility=NULL: the DNA accessibility data. If NULL, it assumes that the entire genome is accessible
### setSequence=NULL: a list of GenomicRaegions where the occupancy will be computed. If NULL, the occupancy is computed genome wide
### lambda=1: the scalling factor between the PWM and the binding energy E=w*(1/lambda).
### boundMolecules=c(1): a vector containing the number of bound molecules to the DNA that
### norm=FALSE: a logical value indicating whether the occupancy data should be normalised.
### outputBoundMoleculesOccupancy=TRUE: logical value indicating whether to compute the occupancy based on PWM score and TF abundance or not.
### backgroundSignal=NULL: the background signal (mean or median of the ChIP-seq data)
### maxSignal=NULL: the highest peak in the ChIP-seq data.
### outputChIPseqProfile=FALSE: logical value indicating whether the method should output the ChIP-seq profile or not.
### chipMean=150: the mean of the ChIP peak see Kaplan et al. 2011
### chipSd=150: the standard deviation of the ChIP peak 
### chipSmooth=NULL: the size of the window used to smooth the landscape
### removeBackground=0: the maximum number of reads to remove
### profile=NULL: ChIP-seq profile
### return: a list containing the following elements:
###			occupancyAbundance: a list with the occupancy based on the PWM score and TF abundance
###			occupancyAbundanceChIP: the same as occupancyAbundance, except that the probability of occupancy display a ChIP-seq shape
###			correlation: a list with the correlation between occupancyAbundanceChIP and the profile
###			MSE: a list with the mean square error between occupancyAbundanceChIP and the profile
###			meanCorrelation: a list with the mean correlation between occupancyAbundanceChIP and the profile for each TF
###			meanMSE: a list with the average mean square error between occupancyAbundanceChIP and the profile for each TF
###			meanTheta: a list with the mean theta value for each TF
computeOccupancy <-function(TF, PWMScore, indexPWMThresholded, DNALength, averageExpPWMScore, DNAAccessibility=NULL, setSequence=NULL, lambda=1, boundMolecules=c(1), norm=FALSE, outputBoundMoleculesOccupancy=TRUE, backgroundSignal=NULL, maxSignal=NULL, outputChIPseqProfile=FALSE, chipMean=150, chipSd=150, chipSmooth = NULL, removeBackground = 0, profile=NULL){
	

	if(is.null(backgroundSignal) | length(TF)!=length(backgroundSignal)){
		backgroundSignal=rep(0,times=length(TF));
	}
	
	
	if(is.null(maxSignal) | length(TF)!=length(maxSignal)){
		maxSignal=rep(0,times=length(TF));
	}	

	occupancyAbundance = vector("list",length(TF));
	occupancyAbundanceChIP = vector("list",length(TF));
	names(occupancyAbundance)=TF;
	names(occupancyAbundanceChIP)=TF;
		
	correlation = NULL;
	MSE=NULL;		
	if(!is.null(profile) && length(profile)== length(TF)){
		correlation = vector("list",length(TF));
		names(correlation)=TF;
		MSE = vector("list",length(TF));
		names(MSE)=TF;
		meanCorrelation = vector("list",length(TF));
		names(meanCorrelation)=TF;
		meanMSE = vector("list",length(TF));
		names(meanMSE)=TF;
		meanTheta = vector("list",length(TF));
		names(meanTheta)=TF;

	}
		
		
		
	#compute occupancy
	for(i in 1:length(TF)){
		occupancyAbundanceLocal = NULL;
		occupancyAbundanceChIPLocal = NULL;
		index=1;
		if(outputBoundMoleculesOccupancy){
			print(paste("compute occupancy ",TF[i],sep=""));
			occupancyAbundanceLocal = computeOccupancyPWMAbundance(PWMScore=PWMScore[[i]], DNALength=DNALength, averageExpPWMScore=averageExpPWMScore[[i]], DNAAccessibility=DNAAccessibility,
										indexPWMThresholded=indexPWMThresholded[[i]], boundMolecules=boundMolecules[i],lambda=lambda[i], norm=norm,
										backgroundSignal=backgroundSignal[i], maxSignal=maxSignal[i]);									
			# output the ChIP-seq like profiles
			if(outputChIPseqProfile){
				occupancyAbundanceChIPLocal=occupancyAbundanceLocal;
				if(!is.null(correlation)){
					correlation[[i]] = rep(0,length(occupancyAbundanceChIPLocal));
					MSE[[i]] = rep(0,length(occupancyAbundanceChIPLocal));
				}
				for(j in 1:length(occupancyAbundanceChIPLocal)){
					occupancyAbundanceChIPLocal[[j]][which(occupancyAbundanceChIPLocal[[j]] < removeBackground)]=0;
					occupancyAbundanceChIPLocal[[j]]=generateChIPProfile(occupancyAbundanceChIPLocal[[j]],chipMean,chipSd,chipSmooth,TRUE);
					if(!is.null(correlation)){
						correlation[[i]][j] = cor(profile[[i]][[j]],occupancyAbundanceChIPLocal[[j]]);
						MSE[[i]][j] = sum((profile[[i]][[j]]-occupancyAbundanceChIPLocal[[j]])^2) / length(profile[[i]][[j]]);
					}	
				}
				
				meanMSE[[i]]=mean(MSE[[i]][which(!is.na(MSE[[i]]))])*1000;
				meanCorrelation[[i]]=mean(correlation[[i]][which(!is.na(correlation[[i]]))])
				meanTheta[[i]]=meanCorrelation[[i]]/meanMSE[[i]];
			}
				
		} 				
		occupancyAbundance[[i]] = occupancyAbundanceLocal
		occupancyAbundanceChIP[[i]] =occupancyAbundanceChIPLocal	
	}
	
	occupancy = list("occupancyAbundance"=occupancyAbundance, "occupancyAbundanceChIP"=occupancyAbundanceChIP, "correlation"=correlation, "MSE"=MSE, "meanMSE"=meanMSE, "meanCorrelation"=meanCorrelation,"meanTheta"=meanTheta);
	
	return(occupancy);
	
}

### Computes the occupancy, given a list of PWM scores and the number of bound molecules. 
### PWMScore: a list of vectors with the PWM scores at each position on the genome. The binding energy (E) can be approximated by the PWM score (w) as E=-w
### DNALength: the length of the DNA
### averageExpPWMScore: the average exponential PWM score
### DNAAccessibility=NULL: a vector with the accessibility level at each position. The accessibility needs to be normalised to its highest value so that it resides win the interval [0,1]. If this is NULL accessibility is not used to compute the occupancy profile. 
### indexPWMThresholded=NULL: the sites that have a PWM score higher than a threshold, if NULL all sites are considered
### boundMolecules=1: the number of bound molecules to the DNA
### lambda=1: the scalling factor between the PWM and the binding energy E=w*(1/lambda).
### norm=TRUE: if this is TRUE the profile is normalised by the highest value. 
### backgroundSignal=0: the background signal (mean or median of the ChIP-seq data)
### maxSignal = 1: the highest peak in the ChIP-seq data.
### return: a list consisting of a vector with the occupancy for each locus
computeOccupancyPWMAbundance <- function(PWMScore, DNALength, averageExpPWMScore, DNAAccessibility=NULL, indexPWMThresholded=NULL, boundMolecules=1, lambda=1, norm=TRUE, backgroundSignal=0, maxSignal=1){	
	result=NULL;
	
	
	result = vector("list",length(PWMScore));
	names(result)=names(PWMScore);
	for(i in 1:length(PWMScore)){
		result[[i]] = rep(0,length(PWMScore[[i]]));
		if(!is.null(DNAAccessibility) & is.list(DNAAccessibility) & length(PWMScore)==length(DNAAccessibility)){
			minLength = min(length(PWMScore[[i]]),length(DNAAccessibility[[i]]));	
			if(is.null(indexPWMThresholded)){
				indexes=1:minLength;
			} else{
				indexes=indexPWMThresholded[[i]][which(indexPWMThresholded[[i]]<minLength)];
			}
			
			
			result[[i]][indexes] = (boundMolecules*DNAAccessibility[[i]][indexes]*exp((1/lambda)*PWMScore[[i]][indexes]))/(boundMolecules*DNAAccessibility[[i]][indexes]*exp((1/lambda)*PWMScore[[i]][indexes])+(((DNALength))/(1))*averageExpPWMScore);
		} else{
			if(is.null(indexPWMThresholded)){
				indexes=1:length(PWMScore[[i]]);
			} else{
				indexes=indexPWMThresholded[[i]];
			}
			result[[i]][indexes] = 1/(1+(((DNALength-(boundMolecules-1)))/(boundMolecules))*averageExpPWMScore*exp(-(1/lambda)*PWMScore[[i]][indexes]));
		}
		
		result[[i]] = backgroundSignal + result[[i]]*(maxSignal-backgroundSignal);
		
	} 
	
#normalise the profile
	if(norm){
		maxOccupancy = max(c(maxSignal,max(unlist(lapply(result,max)))));
		result=lapply(result,"/",maxOccupancy);			
	}	


	return(result)
}


### genarates an array with the names of the loci
### setSequence. A GenomicRanges object with the list of all loci
### setSequenceName: The default names of the loci 
generateLociName <- function(setSequence, setSequenceName){
	if(!is.null(setSequenceName) & length(setSequence)==length(setSequenceName)){
		names=setSequenceName;
	} else{
		names=paste(seqnames(setSequence),":",start(ranges(setSequence)),"..",end(ranges(setSequence)),sep="");
	}
	return(names);
}


### Computes an array with the set of parameters thatfits best the data
### occupancy: a list (lambda) of lists (bound mole) of the dataset generated by computeOccupancy method
### TF: the array of the names of TFs for which the heatmaps will be generated.
### lambdas:an array with all used values of lambda when computing the occupancy lists.
### boundMolecules: an array with all used values of number of bound molecules when computing the occupancy lists.
### parameter="theta": a string indicating which parameter will be used when quantifing the fit. Can be "theta", "correlation", "MSE"
### returns a matrix consisting of two columns (one for lambda and one for the number of bound molecules) a a number of rows eaqual to the number of TFs.
getOptimalSetOfParameters <-function(occupancy, TF, lambdas, boundMolecules,parameter="theta"){
	
	results=matrix(0, nrow=length(TF),ncol=2);
	
	rownames(results)=TF;
	colnames(results)=c("lambda","boundMolecules");
	
	for(TFid in 1:length(TF)){	

		
		if(parameter=="correlation"){
			meanCorrelation=matrix(0,nrow=length(lambdas),ncol=length(boundMolecules));
			for(j in 1:length(lambdas)){
				for(k in 1:length(boundMolecules)){
					meanCorrelation[j,k] = occupancy[[j]][[k]]$meanCorrelation[[TF[TFid]]];
				}
			}
			optimalIndex=which(meanCorrelation==max(meanCorrelation), arr.ind = TRUE);			
		} else if(parameter=="MSE"){
			meanMSE=matrix(0,nrow=length(lambdas),ncol=length(boundMolecules));
			for(j in 1:length(lambdas)){
				for(k in 1:length(boundMolecules)){
					meanMSE[j,k] = occupancy[[j]][[k]]$meanMSE[[TF[TFid]]];
				}
			}
			optimalIndex=which(meanMSE==min(meanMSE), arr.ind = TRUE);			
		} else{
			meanTheta=matrix(0,nrow=length(lambdas),ncol=length(boundMolecules));
			for(j in 1:length(lambdas)){
				for(k in 1:length(boundMolecules)){
					meanTheta[j,k] = occupancy[[j]][[k]]$meanTheta[[TF[TFid]]];
				}
			}
			optimalIndex=which(meanTheta==max(meanTheta), arr.ind = TRUE);
		}
		results[TFid,1]=optimalIndex[1];
		results[TFid,2]=optimalIndex[2];
	}
	
	return(results);
}


