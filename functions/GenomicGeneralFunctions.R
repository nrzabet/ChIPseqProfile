
### process the read density from DNaseI data and compute the probability that a region is accessible. see Kaplan et al. 2011
### readDensity: the number of reads at each postion
### alpha a factor to compute probability
### beta a factor to compute probability
processReadDensity <- function(readDensity, alpha, beta){
	prob = 1/(1+exp(-beta*readDensity + alpha));
	return(prob);
}

### extract the set of genes from a tsv file
### filename: the file storing the gene reference
### return: a list containing 4 GRanges objects. The GRanges objects contain the coordinates of the exon, intron 5'UTR and 3'UTR.
extractGeneList <- function(filename){
	buffer=read.table(filename,stringsAsFactors=FALSE);
	
	indexExon=which(buffer[,5]=="exon");
	indexIntron=which(buffer[,5]=="intron");
	index5UTR=which(buffer[,5]=="five_prime_UTR");
	index3UTR=which(buffer[,5]=="three_prime_UTR");
	
	grExon = GRanges(seqnames=Rle(buffer[indexExon,1]), ranges=IRanges(buffer[indexExon,2],buffer[indexExon,3]), strand=buffer[indexExon,4]);
	names(grExon)=buffer[indexExon,6];
	
	grIntron = GRanges(seqnames=Rle(buffer[indexIntron,1]), ranges=IRanges(buffer[indexIntron,2],buffer[indexIntron,3]), strand=buffer[indexIntron,4]);
	names(grIntron)=buffer[indexIntron,6];
	
	gr5UTR = GRanges(seqnames=Rle(buffer[index5UTR,1]), ranges=IRanges(buffer[index5UTR,2],buffer[index5UTR,3]), strand=buffer[index5UTR,4]);
	names(gr5UTR)=buffer[index5UTR,6];

	gr3UTR = GRanges(seqnames=Rle(buffer[index3UTR,1]), ranges=IRanges(buffer[index3UTR,2],buffer[index3UTR,3]), strand=buffer[index3UTR,4]);
	names(gr3UTR)=buffer[index3UTR,6];

	gr = list("exon"=grExon, "intron"=grIntron, "5UTR"=gr5UTR, "3UTR"=gr3UTR);
	gr
}





### extract the set of enhancers locations from a bed file
### filename: the file storing the enhancer list
### return: a GRange object
extractEnhancers <- function(filename){
	buffer=read.table(filename,stringsAsFactors=FALSE);
	
	enhancer = GRanges(seqnames=Rle(buffer[,1]), ranges=IRanges(buffer[,2],buffer[,3]), strand=rep('*',times=nrow(buffer)));
	names(enhancer)=buffer[,4];
	return(enhancer)
}


### extract a list with the occupancy at each loci from a ChIP-seq data set
### profile: the matrix with the ChIP-seq data with the following columns: chromosome, start, end, value
### setSequence: a GenomicRange object with all loci where the ChIP-seq data is extracted at
### maxSignal=NULL:the value to which the data is normalised. If NULL no normalisation is performed 
### removeBackground=0: The threshold under which all data is assume to be 0. Can be used to truncate noisy data
### chipSmooth=NULL: The window over which to smooth the data. If NULL no smooting is performed
### return: a list consistsing of an array at each locus with the normalised ChIP-seq signal
extractOccupancyDataAtLoci <- function(profile, setSequence, maxSignal=NULL, removeBackground=0, chipSmooth=NULL){
	#compute overlayed occupancy
	
	# create a list with the occupancy at each loci 
	occupancyNorm = vector("list",length(setSequence));
	if(!is.null(setSequence) & length(names(setSequence))>0){
		names(occupancyNorm)=names(setSequence);
	} else{
		names(occupancyNorm)=paste(seqnames(setSequence),":",start(ranges(setSequence)),"..",end(ranges(setSequence)),sep="");
	}
	
	
	for(gr in 1:length(setSequence)){		
		chr=as.vector(seqnames(setSequence))[gr];
		range=ranges(setSequence);
		startPositon=as.vector(range@start)[gr];
		endPosition = as.vector(end(range))[gr];
		startPos=startPositon;
		endPos=endPosition;
		positions=startPos:endPos; 
		occupancyNorm[[gr]]=rep(0,(endPos-startPos+1));
		
		if(!is.null(profile) & length(profile)>0){
			profile[,2] = as.numeric(profile[,2]);
			profile[,3] = as.numeric(profile[,3]);
			profile[,4] = as.numeric(profile[,4]);
			indexes = intersect(intersect(which(profile[,2]<endPos), which(profile[,3]>=startPos)), which(profile[,1]==chr));
			if(length(indexes)>0){
				
				#extract the occupancy
				for(sector in indexes){
					minPos = max(profile[sector,2]-startPos,1);
					maxPos = min(profile[sector,3]-startPos,length(occupancyNorm[[gr]]));					
					occupancyNorm[[gr]][minPos:maxPos] = profile[sector,4];
				}
				
				#normalise the occupancy signal	
				if(is.null(maxSignal)){
					occupancyNorm[[gr]]=occupancyNorm[[gr]]/max(profile[,4]);
				} else{
					occupancyNorm[[gr]]=occupancyNorm[[gr]]/maxSignal;
				}
				
				#remove the occupancy lower than a threshold
				occupancyNorm[[gr]][which(occupancyNorm[[gr]] < removeBackground)]=0;
				
				# smooth the occupancy vector is required
				if(!is.null(chipSmooth)){
					occupancyNorm[[gr]]=vectorSmooth(occupancyNorm[[gr]],chipSmooth,TRUE);
				}
			}
		}
	}
	return(occupancyNorm);
}


### Smooths a vector
### x: the vector
### smooth=NULL: the window size used to smooth the profile, if NULL then the profile is not smoothed
### norm=FALSE: this is true if the vector is normized at the end to the maximum value 
### return: a smoothed vector
vectorSmooth <- function(x, smooth=NULL,norm=FALSE) {
	originalMax=max(x);
	result=x;
	if(!is.null(smooth)){
        if((smooth %% 2) == 0){smooth = smooth - 1}
        mid = round(smooth/2,0) + 1
        d = smooth - mid
        for(i in mid:(length(x) - d)) {
            result[i] = mean(x[max(0,(i-d)):min(length(x),(i+d))]);
        }
    }
	
	
	if(norm){
		result=originalMax*result/max(result);
	}	

	result
}


### Computes the nucleotide frequency in DNA sequence
### sequences: the DNAStringSet object with the list of DNA sequences
### return: a 4 elements vector of the nucleotide frequency, with rows: A, C, G, T
computeBPFrequency <- function(DNAsequence){
	
	BPFreq=apply((alphabetFrequency(DNAsequence, baseOnly=TRUE)),2,sum)[1:4]
	
	#normalise
	BPFreq=BPFreq/sum(BPFreq);
	
	return(BPFreq);
}


### Computes the PWM matrix from the PFM matrix
### PFM: the PFM matrix from file
### pseudocount=1: the pseudocount used to compute the PWM
### BPFrequency=rep(0.25,4): a vector with the genome wide frequencies of the four nucleotides.
### naturalLog=FALSE: if this is TRUE the method uses the natural logarithm  when computing the PWM (Berg and von Hippel, 1987), otherwise the method uses logarithm in base 2 (Stormo, 2000)
### noOfSites=NULL: the number of DNA sequences used to compute the PFM. If NULL it will be computed based on the PFM
### return: a 4 row matrix of the PWM, with rows: A, C, G, T
computePWM <- function(PFM=NULL, pseudocount=1, BPFrequency=rep(0.25,4),naturalLog=FALSE,noOfSites=NULL){
	if(min(PFM)<0){
		#the provided PFM is the actual PWM
		PWM=PFM;
	} else{
		#get the sum on each colum (it should be equal)
		if(is.null(noOfSites)){
			noOfSites = max(apply(PFM,2,sum));
		}

		#normalise the pfm
		for(i in 1:ncol(PFM)){
			PFM[,i]= PFM[,i]/max(apply(PFM,2,sum));
		}
		
		#generate the pwm based on the pfm
		PWM = computePWMFromPFM(PFM,max(noOfSites), pseudocount, BPFrequency,naturalLog);	
	}
	
	
	return(PWM);
}

### Computed the PWM matrix using the PFM matrix
### PFM: the PFM matrix
### noOfSites: the total number of sites used to compute the PFM matrix
### pseudocount=1: the pseudocount used to compute the PWM
### BPFrequency=rep(0.25,4): a vector with the genome wide frequencies of the four nucleotides .
### naturalLog=FALSE: if this is TRUE the method uses the natural logarithm  when computing the PWM (Berg and von Hippel, 1987), otherwise logarithm  in base 2 is used (Stormo, 2000)
### return: a 4 row matrix of the PWM, with rows: A, C, G, T
computePWMFromPFM <- function(PFM,noOfSites, pseudocount=1, BPFrequency=rep(0.25,4),naturalLog=FALSE){
	PWM = matrix(0,nrow=nrow(PFM),ncol=nrow(PFM));
	PWM = PFM;
	for(i in 1:nrow(PFM)){
		for(j in 1:ncol(PFM)){
			PWM[i,j]=(PFM[i,j]*noOfSites+BPFrequency[i]*pseudocount)/(noOfSites+pseudocount);
			if(naturalLog){
				PWM[i,j]=log(PWM[i,j]/BPFrequency[i]);
			} else{
				PWM[i,j]=log2(PWM[i,j]/BPFrequency[i]);
			}
		}
	}
	return(PWM);
}


### Attempts to parse a file and extract the PFM. It uses the following formats: RAW, TRANSFAC, Jaspar or constructs one from sequences. For the correct formats see http://www.benoslab.pitt.edu/stamp/help.html#input
### filename: the file storing the motif
### format="raw": the format of the PWM: raw, transfac, jaspar, seqs (list of binding sequences)
### return: a PFM matrix, with rows: A, C, G, T
parsePFM <- function(filename, format="raw"){
	PFM=NULL;
	if(format=="raw"){
		PFM=parseRawPFM(filename);
	}
	if(format=="transfac"){
		PFM=parseTransfacPFM(filename);
	}
	if(format=="JASPAR"){
		PFM=parseJASPARPFM(filename);
	}
	if(format=="sequences"){
		PFM=parsePFMFromFile(filename);
	}
	
	#transpose the PFM/PWM if required
	if(nrow(PFM)!=4){
		if(ncol(PFM) == 4){
			PFM=t(PFM);
		} else{
			stop("Could not load PFM/PWM of 4 rows or columns. Wrong PFM/PWM");
		}
	}
	
	
	return(PFM);
}


### Extracts the PFM matrix from a raw PSSM style file. It also works for MEME output
### filename: the file containing the motif
### return: the PFM matrix
parseRawPFM <- function(filename){
	buffer = read.table(filename,row.names=NULL,stringsAsFactors=FALSE);
	if(ncol(buffer)==4){
		PFM=t(buffer);
		rownames(PFM)=c("A","C","G","T");
	} else{
		PFM=NULL;
	}
	
	return(PFM);	
}

### Extracts the PFM matrix from a transfac style file.
### filename: the file containing the motif
### return: a 4 row PFM matrix, with rows: A, C, G, T
parseTransfacPFM <- function(filename){
	buffer = read.table(filename,row.names=NULL,stringsAsFactors=FALSE);
	if(ncol(buffer)>4){
		PFM=t(buffer[,2:5]);
		if(is.character(PFM[,1])){
			PFM=pfm[,2:ncol(PFM)];
		}
		rownames(PFM)=c("A","C","G","T");
	} else{
		PFM=NULL;	
	}
	return(PFM);	
}

### Extracts the PFM matrix from a jaspar style file
### filename: the file containing the motif
### return: a 4 row PFM matrix, with rows: A, C, G, T
parseJASPARPFM <- function(filename){
	buffer = readLines(filename);
	error=FALSE;
	for(i in 1:length(buffer)){
		if(!error){
			buffer[i]=sub("[","",buffer[i],fixed=TRUE);
			buffer[i]=sub("]","",buffer[i],fixed=TRUE);
			cells = strsplit(buffer[i]," +",fixed=FALSE)[[1]];
			
			if(length(cells)>1){
				values=as.numeric(cells[2:(length(cells))]);
				if(i==1){
					PFM=values;
				} else{
					PFM=rbind(PFM,values)
				}
			} else{
				error=TRUE;
			}
		}
	} 
	
	if(error){
		PFM=NULL;
	} else{
		rownames(PFM)=c("A","C","G","T");
	}
	return(PFM);	
}




### Reads a list of sequences from a file and produces the PFM mtrix;
### fileame: the filename with the sequences
### return: a 4 row PFM matrix, with rows: A, C, G, T
parsePFMFromFile <- function(fileame){
	sequences=readDNAStringSet(fileame);	
	if(!is.null(sequences) & length(sequences)>0){
		PFMExtended = consensusMatrix(sequences,as.prob=TRUE, baseOnly=TRUE);
		PFM=PFMExtended[1:4,];
	} else{
		PFM=NULL;
	}
	
	return(PFM);
}





### Computes the PWM score of a DNAStringSet given a PWM matrix.  
### PWM: the PWM matrix consisting of 4 rows corresponding to A, C, G, T
### DNASequence: the DNAStringSet object with the list of DNA sequences
### strand="+": the strand used to compute the PWM score; can be "+", "-", "+-"
### strandRule="max": the rule used to combine the PWM score on two strands; can be "max", "sum" or "mean"
### return: a list consisting of a vector with the PWM scores for each DNAString
scoreDNAStringSet <-function(PWM, DNASequenceSet, strand="+", strandRule="max"){
	scoreSet = vector("list", length(DNASequenceSet));
	names(scoreSet) = names(DNASequenceSet);
	for(i in 1:length(DNASequenceSet)){
		sequenceLength = (length(DNASequenceSet[[i]])-ncol(PWM)+1);
		scorePositive = NULL;
		scoreNegative = NULL;
		if(length(grep("+",strand))){
			scorePositive = PWMscoreStartingAt(PWM,DNASequenceSet[[i]],starting.at=1:sequenceLength);
		}
		
		if(length(grep("-",strand))){
			scoreNegative = rev(PWMscoreStartingAt(PWM,reverseComplement(DNASequenceSet[[i]]),starting.at=1:sequenceLength));
		}
					
		if(is.null(scorePositive) & is.null(scoreNegative)){
			# no strand is specified so consider the positive strand
			scoreSet[[i]] = PWMscoreStartingAt(PWM,DNASequenceSet[[i]],starting.at=1:sequenceLength);
		} else if(!is.null(scorePositive) & !is.null(scoreNegative)){
			# both strands are specified 
			if(strandRule =="sum"){
				scoreSet[[i]] = scorePositive + scoreNegative;
			} else if(strandRule =="mean"){
				scoreSet[[i]] = (scorePositive + scoreNegative)/2;
			} else{
				scoreSet[[i]] = pmax(scorePositive,scoreNegative);
			} 
		} else{
			if(!is.null(scoreNegative)){
				scoreSet[[i]] = scoreNegative;
			} else{
				scoreSet[[i]] = scorePositive;
			}
		}
	}
	return(scoreSet);
}



### Generates an in silico ChIP-seq profile
### inputVector: the occupancy vector at nucleotide resolution
### mean: the mean of the width the ChIP-seq peak 
### sd: the sd of the width of the ChIP-seq peak
### smooth: the window size used to smooth the profile, if NULL then the profile is not smoothed
### norm: if this is TRUE the profile is normalised by the highest value.
### quick: if this is TRUE the ChIP peak is computed with an approximation  
### return: a vector with the value at each position representing the in-silico ChIP-seq signal.
generateChIPProfile <- function(inputVector, mean, sd, smooth = NULL, norm = TRUE, quick = TRUE) {
	originalMax=max(inputVector);
	
	var = sd^2
	shp = mean^2/var
	scl = var/mean
	l = length(inputVector)

	peakSignificantSize=1250;
	
	f = dgamma(0:length(inputVector), shape = shp, scale = scl)
	F = rev(cumsum(rev(f)))
	sp = c(rev(F[2:peakSignificantSize]), F[1:peakSignificantSize])
	lp = length(sp)
	hp = (lp - 1)/2

	peak.centres = which(inputVector > mean(inputVector))
	peaks = vector("numeric", l)
	
	if(!quick){	
		for(pc in peak.centres) {
			thisPeak = vector("numeric", l)
			thisPeak[pc:l] = F[1:(l-pc+1)]
			if(pc!=1) thisPeak[1:(pc-1)] = F[pc:2]
			peaks = peaks + thisPeak * inputVector[pc]
		}
	} else {
		for(pc in peak.centres) {
			thisPeak = sp
			left = pc - hp
			right = pc + hp
			cutLeft = max(0, 0 - left)
			cutRight = -min(0, l - right)	
			if(cutLeft != 0){
				thisPeak = thisPeak[-(1:(cutLeft+1))]
			}
			if(cutRight != 0){
				thisPeak = thisPeak[-((lp-cutRight+1):lp)]
			}
			left = max(1, left)
			right = left + length(thisPeak) - 1
			peaks[left:right] = peaks[left:right] + thisPeak * inputVector[pc]
		}	
	}

	if(!is.null(smooth)){
		if((smooth %% 2) == 0){smooth = smooth - 1}
		mid = round(smooth/2,0) + 1
		d = smooth - mid
		for(i in mid:(length(peaks) - d)) {
			peaks[i] = mean(peaks[(i-d):(i+d)])
		}
	}

	
	if(norm && max(peaks)!=0){
		peaks=originalMax*peaks/max(peaks);
	}	
	
	return(peaks)
}



