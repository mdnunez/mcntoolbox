## svderp.R provides an example of how to use singular value decomposition (SVD) to find single-trial event-related potentials (ERPs)

##Embedded in the script is a simple function to compute single-trial ERPs in R

## Citation: 
# Nunez M.D., Vandekerchove, J., Srinivasan, R. (2016) How attention influences perceptual decision making: Single-trial EEG correlates of drift-diffusion model parameters. Journal of Mathematical Psychology.

## Copyright 2015 Michael D. Nunez

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
 #   along with this program.  If not, see <http://www.gnu.org/licenses/>.

###Record of Revisions###
#   Date           Authors                   Description of change
#   ====        =================            =====================
#  12/11/15      Michael D. Nunez                Original Code

#######
###This function calculates single trial ERPs using the SVD method###
svderp <- function(eeg,trialmarks,trialsamps,baseline=NULL,comp=1) {
#eeg is a vector that marks the point of reference of each trial
#trialsamps + trialmarks[i] gives the indices of specific trials 
#baseline + trialmarks[i] are the indices for baselining each trial
#comp is the component number to compute, comp = 1 will have the highest percent of variance explained

ntrials <- length(trialmarks)
eegsum <- array(data=0,dim=c(length(trialsamps),dim(eeg)[2]))

#Calculate ERP
if(is.null(baseline)) {
	#No baseline
	for (i in 1:ntrials) {
	    eegsum <- eegsum + eeg[trialmarks[i] + trialsamps,]
	}
	} else {
	#Baseline
	for (i in 1:ntrials) {
	    baseline <- t(as.array(apply(eeg[trialmarks[i] + trialsamps,],2,mean)))
	    eegsum <- eegsum + eeg[trialmarks[i] + trialsamps,] - baseline[rep(1,length(trialsamps)),]
	}
}

erp <- eegsum/ntrials #The ERP is the mean of EEG across trials, time-locked to specific events

#SVD of ERP (equivalent to PCA)
svdouts <- svd(erp)

#Percent of variance explained by each principal component
perexp = (svdouts$d^2)/sum(svdouts$d^2);

#nchan*1 weights for each electrode is just the right singular vector for the specific component
#Note that the weight vector's sign is arbitrary, it may have to be changed to reflect the ERP
weights <- svdouts$v[,comp]


#Single-trial ERPs
sterp = eeg %*% weights

#Use trialmarks to find single-trial ERPs in svderp$sterp via svderp$sterp[trialmarks(i) + trialsamps,]

#Return ERP, percent of variance explained by component comp, nchan*1 weight vector, and single-trial erps
outs <- list("erp" = erp, "pexp" = perexp, "w"=weights,"sterp" = sterp)
return(outs)

}
#######

##The following code simulates EEG data with an embedded response to external stimuli and compares the raw EEG at one channel with a high signal-to-noise ratio to the single-trial ERP

set.seed(380)

ntrials <- 100
nchans <- 128
triallen <- 1000

trialt <- 1:triallen

#Simulated ERP
trueerp <- 1000*sin(trialt*2*pi*(4/triallen))*dnorm(trialt,mean=triallen/5,sd=triallen*.15)-sin(trialt*2*pi*(1/triallen))

#Simulated EEG
eeg <- array(dim = c(triallen*ntrials,nchans))

chanweights <- rep(0,1,nchans)
chanweights[21:30] <- -c(rep(5,5),rep(10,5))
chanweights[81:90] <- c(rep(5,5),rep(10,5))

for (k in 1:ntrials) {
	for (c in 1:nchans) {
		eeg[(1:triallen)+(k-1)*triallen,c] <- trueerp*chanweights[c] + rnorm(triallen,sd = 30)
	}
}

trialmarks <- seq(1,triallen*ntrials,triallen) #This vector marks the point of reference of each trial
trialsamps <- 0:(triallen-1) #This vector + trialmarks[i] gives the indices of specific trials 

svdouts <- svderp(eeg,trialmarks,trialsamps,comp=1) #Compute the single-trial ERP


####Plots
par(mfrow=c(3,1))
#Plots the first trial of EEG of one of the best electrodes for the stimulus response
plot(1:triallen,eeg[1:triallen,90],'l',xlab='',ylab=expression(paste(mu,"V")),main = '1st Trial and Best Channel of Raw EEG',cex.lab=2,cex.axis=2,cex.main=2,cex.sub=2)
#Plots the ERP at one of the best electrodes
plot(1:triallen,svdouts$erp[,90],'l',xlab='',ylab=expression(paste(mu,"V")),main = 'Best Channel of ERP',cex.lab=2,cex.axis=2,cex.main=2,cex.sub=2)
#Plots the first trial of the single-trial ERP
plot(1:triallen,svdouts$sterp[1:triallen,],'l',xlab='Time Post-Stimulus (ms)',ylab=expression(paste(mu,"V")),main = '1st Trial of Single-trial ERP',cex.lab=2,cex.axis=2,cex.main=2,cex.sub=2)

