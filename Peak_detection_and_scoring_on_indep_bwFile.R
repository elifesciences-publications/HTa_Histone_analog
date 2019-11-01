#This script is dedicated to the detection of peaks from one bigwig file, and their scoring according another bigwig file (MN1 and MN2)

rm(list=ls())

library(rtracklayer);library(ggplot2)
#This script use a modified version of the NucleR package https://bioconductor.org/packages/release/bioc/html/nucleR.html
source("/Users/ahocher/Documents/Scripts/R/nucleR_2.12.1_AH_edited_asy/R/filterFFT.R")
source("/Users/ahocher/Documents/Scripts/R/nucleR_2.12.1_AH_edited_asy/R/helpers.R")
source("/Users/ahocher/Documents/Scripts/R/nucleR_2.12.1_AH_edited_asy/R/peakDetection.R")
source("/Users/ahocher/Documents/Scripts/R/nucleR_2.12.1_AH_edited_asy/R/peakScoring.R")







#Setting the fourier transform filtering value for all downstream analysis
PcVal=0.012
Zthr=0.25

#Here we set the length on which the second mnase track will be scored :
DyadLength=80

#Loading of MNase data : 
MN1file="/Users/ahocher/Documents/TEAC_HTa_Study/Results/MNase_Experiments/Sequencing/2018_07_PE_TEAC/Coverage_tracks_Merged_reads_Bowtie2/RGC_Input_Zscore_Normalized/Merged_replicates/Merged_replicates_day2_40_65bp.bw_Zscore_normalized.bw"
MN1=import.bw(MNShortfile,as="RleList")

MN2file="/Users/ahocher/Documents/TEAC_HTa_Study/Results/MNase_Experiments/Sequencing/2018_07_PE_TEAC/Coverage_tracks_Merged_reads_Bowtie2/RGC_Input_Zscore_Normalized/Merged_replicates/Merged_replicates_day2_70_100bp.bw_Zscore_normalized.bw"
MN2=import.bw(MNSmallfile,as="RleList")


#Because MNase data is under the form of Z-score we have to add a value to all values to make it positive
MNCov=MN1
MNCovdt <- as.data.frame(MNCov)
MinVal=min(MNCovdt$value)
MNCov=MNCov-MinVal

FilterMN  <-  filterFFT(data=MNCov, pcKeepComp=PcVal, showPowerSpec=F, useOptim=T, mc.cores=4)

#Checking which value corresponds to z-score = 0.25, necessary as the fourier filtering induced changes in the values. 
Thr=Zthr*sd(FilterMN$NC_002578)+mean(FilterMN$NC_002578)


#In opposition to previous script we first detect peak without scoring them.
LesPeaksCoord <- peakDetection(FilterMN, threshold=Thr, score=F,width=90,dyad.length = DyadLength,min.cov = Thr)


#Transforming short fragments bw for peakscoring based on those ( the format has to be a vector for some reason)
MNCovVect=vector(mode = "list")
MNCovVect$NC_002578=as.data.frame(MN2)$value
LesPeaks=peakScoring(peaks = LesPeaksCoord ,data = MNCovVect ,threshold = Thr, dyad.length = DyadLength)


LesPeaksdt <- as.data.frame(LesPeaks)

LesPeaksdt=LesPeaksdt[order(LesPeaksdt$score_asy),]


#Export a specific file as bed for further use on DNA sequence plot :
ToexportBed=LesPeaksdt[,c(1,2,3,8,9)]
names(ToexportBed)[c(4,5)]=c("name","score")
export.bed(ToexportBed,"Peaks_detected_on_MN1_Scored_on_MN2.bed")


#Annotation and export of separate bed files based on scores :


Peaks=import.bed("Peaks_detected_on_MN1_Scored_on_MN2.bed")

Peaks=Peaks[order(start(Peaks))]
strand(Peaks[Peaks$score>1])="-"
strand(Peaks[Peaks$score<1])="+"

IntensePeaks=Peaks[as.numeric(Peaks$name)>0.7]
FaintPeaks=Peaks[as.numeric(Peaks$name)<0.7]

SymmetricalPeaks=IntensePeaks[abs(log2(IntensePeaks$score))<0.25]
ASymmetricalPeaks=IntensePeaks[abs(log2(IntensePeaks$score))>0.25]


export.bed(Peaks,"Stranded_peaks.bed")
export.bed(IntensePeaks,"Intense_peaks_sup07.bed")
export.bed(FaintPeaks,"Faint_peaks_sup07.bed")
export.bed(SymmetricalPeaks,"Symmetrical_Peaks_sup07.bed")
export.bed(ASymmetricalPeaks,"ASymmetrical_Peaks.bed")
