#! R

##AL comments: example call to this script: Rscript ./Scripts/activity_score_highfreq_signal_range.R HD2_FC15884027.markDuplicates/TranscriptionFactors/GTRD_ChIP_Only_1000sites/Androgen.Top1000sites.tss Androgen.Top1000sites.tss HD2_FC15884027.markDuplicates/AccessibilityOutput/HD2_FC15884027.markDuplicates_Accessibility1KSites.txt
##to run this from R studio. First needed to move one of the results folders into the scripts directory (symlink didn't work).
##to run add these commands: 
#setwd('/Volumes/ha_g/user/adoebley/ha_lab_scripts/bin/Ulz_TranscriptionFactorProfiling')
#args <- c('HD2_FC15884027.markDuplicates/TranscriptionFactors/GTRD_ChIP_Only_1000sites/Androgen.Top1000sites.tss', 'Androgen.Top1000sites.tss', 'HD2_FC15884027.markDuplicates/AccessibilityOutput/HD2_FC15884027.markDuplicates_Accessibility1KSites.txt')
##also had to comment out the command line args below
####end AL comments


suppressMessages(library(signal)) #this both loads signal and supresses messages

getLowSignal <- function(coverage_signal) {
   normalized <- coverage_signal / mean(coverage_signal)
   low<-sgolayfilt(normalized,3,1001)
   return(low)
}
getHighSignal <- function(coverage_signal,low_signal) {
   normalized <- coverage_signal / mean(coverage_signal)
   high<-sgolayfilt(normalized,3,51)
   high_adjusted<-(high/low_signal)
   return(high_adjusted)
}
find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE))) # AL comment: takes the difference between adjacent points, then the sign of this difference, then the difference between these signs
    pks <- sapply(which(shape < 0), FUN = function(i){ #AL comment: which(shape <0) gets a list of which indexes in shape are <0
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1) #AL: if z> 0 return z, otherwise return 1
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0)) #checks if the value at position is larger than the values within the window from z to w
    })
     pks <- unlist(pks)
     pks
}
#AL summary of find_peaks: This function is finding all the local maxima and minima (saved in shape).
#it then creates a window around each local maxima (i-m to i+m) 
#if the window runs up against the edge, the window is truncated at the edge.
#It then checks the window. If the local maxima i is the largest value in the window, it is considered a peak.
#this data is never actually used


#AL: need to comment out this line to run it from RStudio:
args<-commandArgs(trailingOnly=TRUE) 
data<-read.csv(args[1],header=TRUE,sep="\t")
sample_name<-args[2]
outfile<-args[3]
data$low<-getLowSignal(data$Mean.Cov)
data$high<-getHighSignal(data$Mean.Cov,data$low)

n<-data$TSS.analyzed[1] #get the number of sites analyzed from the data
range<-max(data$high) - min(data$high)
peaks<-find_peaks(data$high,m=20)
peak_positions = data$Position[peaks]
peak_distance = c(diff(peak_positions))
mean_peak_distance = mean(peak_distance)
median_peak_distance = median(peak_distance)

write(paste(sample_name,n,range,length(peaks),mean_peak_distance,median_peak_distance,sep="\t"),file=outfile,append=TRUE)



