#Building up the Input files for LASSO general linear model modeling of MNase data on MATLAB:
rm(list=ls())

library(rtracklayer)
#This script requires the user to create a BSgenome : 
library(BSgenome.Tacidophilum.NCBI.v1)
library(rtracklayer)
library(dplyr)
library(R.matlab)


#Create working directory if it doesn't exist
WorkingDir="/Users/ahocher/Documents/TEAC_HTa_Study/Results/Modeling_HTa_binding/"

if (!dir.exists(WorkingDir)) dir.create(WorkingDir, recursive = TRUE)
setwd(WorkingDir)

# Linear modeling of MNase signal using DNA sequence in H.volcaniss
library(seqTools)


#Determine which Chromosome to study : 
ChromosomeTacid <- as.character(Tacidophilum$NC_002578)

#Function to compute kmers genome wide (flanks are just here to be able to compute number for all base pairs (including the first and the last):
ComputeKmer <- function(Chromosome,Kmerlength,WindowSize)
{
  TotalLength <- nchar(Chromosome)
  FlanksToAdd=100
  Chromosome_FlankFirst <- substr(Chromosome,1,100)
  Chromosome_FlankLast <- substr(Chromosome,TotalLength-FlanksToAdd,TotalLength)
  
  ChromosomeFlanked=paste(Chromosome_FlankLast,Chromosome,Chromosome_FlankFirst,sep="")
  
  
  
  
  #Where to start an stop the calculations to properly include the flanking sequences that had been added.  : 
  TotalLengthFlank <- nchar(ChromosomeFlanked)
  
  WtS <- 1+FlanksToAdd-(WindowSize-1)/2
  #Here the +100 is to account for the addition of 100 basepairs at first.
  WtE <- FlanksToAdd+TotalLength-(WindowSize-1)/2
  
  #TotalLengthTocompute <- TotalLength - WindowSize - Kmerlength + 2
  
  # Run the calculation (no loop required). "start" argument can take a vector of start positions!
  kmerCount_raw <- countDnaKmers(ChromosomeFlanked, Kmerlength, start = WtS:WtE, width = WindowSize)
  kmerCount_df <- as.data.frame(t(kmerCount_raw))
  
  # Keeping the name of the first bast pair
  kmerCount_names <- paste0(1:TotalLength, "bp")
  rownames(kmerCount_df) <- kmerCount_names
  
  #dim(kmerCount_df);head(kmerCount_df)
  # File export
  write.table(kmerCount_df, file = paste0("Kmer_counting_", Kmerlength, "_letters_Tacidophilum_genome_windowSize", WindowSize, ".txt"), row.names = TRUE, quote = FALSE, sep = "\t")
  
}


#Note that window size has to be an odd number.

ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 1,WindowSize = 21)
ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 2,WindowSize = 21)
ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 3,WindowSize = 21)
ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 4,WindowSize = 21)


ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 1,WindowSize = 51)
ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 2,WindowSize = 51)
ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 3,WindowSize = 51)
ComputeKmer(Chromosome = ChromosomeTacid, Kmerlength = 4,WindowSize = 51)


#Loading the Kmer files :
MonoNucFile21=file.path("Kmer_counting_1_letters_Tacid_genome_windowSize21.txt")
DiNucFile21=file.path("Kmer_counting_2_letters_Tacid_genome_windowSize21.txt")
TriNucFile21=file.path("Kmer_counting_3_letters_Tacid_genome_windowSize21.txt")
QuadriNucFile21=file.path("Kmer_counting_4_letters_Tacid_genome_windowSize21.txt")

MonoNucFile=file.path("Kmer_counting_1_letters_Tacid_genome_windowSize51.txt")
DiNucFile=file.path("Kmer_counting_2_letters_Tacid_genome_windowSize51.txt")
TriNucFile=file.path("Kmer_counting_3_letters_Tacid_genome_windowSize51.txt")
QuadriNucFile=file.path("Kmer_counting_4_letters_Tacid_genome_windowSize51.txt")


#In this case we use two different window sizes but it's not necessary : 
Mo21=read.table(MonoNucFile21,header=T,sep="\t");names(Mo21)=paste(names(Mo21),"_21bp",sep="")
Di21=read.table(DiNucFile21,header=T,sep="\t");names(Di21)=paste(names(Di21),"_21bp",sep="")
Tri21=read.table(TriNucFile21,header=T,sep="\t");names(Tri21)=paste(names(Tri21),"_21bp",sep="")
Quadri21=read.table(QuadriNucFile21,header=T,sep="\t");names(Quadri21)=paste(names(Quadri21),"_21bp",sep="")

Mo=read.table(MonoNucFile,header=T,sep="\t")
Di=read.table(DiNucFile,header=T,sep="\t")
Tri=read.table(TriNucFile,header=T,sep="\t")
Quadri=read.table(QuadriNucFile,header=T,sep="\t")




#####LOADING MNASE DATA
#Input data should be normalized (z-score normalisation is generally ok)
isInputDataNormalized=1



MNFile=file.path("/Users/ahocher/Documents/TEAC_HTa_Study/Results/MNase_Experiments/Sequencing/2018_07_PE_TEAC/Coverage_tracks_Merged_reads_Bowtie2/RGC_Input_Zscore_Normalized/40_65/Merged_Trimmed_r2_d2_t15_GTGAAA_L001_R_001_srt_RPGC_45_60bpNormalizedNormalized_by_Input_top_decileInput_normalized_Zscore.bw")
if (isInputDataNormalized==0){
  
  MNdata=import.bw(MNFile,as="RleList")
  
  #Produce a normalized MNase track as the one used for modelling for latter use:
  MNasedt=as.data.frame(MNdata$NC_002578)
  MNase_mean <- mean(MNasedt$value)
  MNase_sd <-sd(MNasedt$value)
  MNase_NormZ=(MNdata-MNase_mean)/MNase_sd
  MNase_NormZ_out <- file.path( paste("/Users/ahocher/Documents/TEAC_HTa_Study/Results/MNase_Experiments/Sequencing/2018_07_PE_TEAC/Coverage_tracks_Merged_reads_Bowtie2/RGC_Input_Zscore_Normalized/",sub(".bw", "", basename(MNFile)),"Input_normalized_Zscore",".bw",sep=""))
  export(MNase_NormZ ,MNase_NormZ_out,format="bw")
  
  MNFile=paste(sub(".bw", "", basename(MNFile)),"Input_normalized_Zscore",".bw",sep="")
  rm(MNdata)
  #Erase and re-import the data out of the loop
  
}
MNdata=import.bw(MNFile,as="RleList")
MN=as.data.frame(MNdata$NC_002578)




# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --
#Choose which part of the genome to use as training and validation sets :
# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --

#Training set : 
Selecta=c(1:round(length(MNdata$NC_002578)/3))

#Test set : 
ASelecta=c((round(length(MNdata$NC_002578)/3)+1):length(MNdata$NC_002578))

InputTrain=cbind(Mo[Selecta,],Di[Selecta,],Tri[Selecta,],Quadri[Selecta,],Mo21[Selecta,],Di21[Selecta,],Tri21[Selecta,],Quadri21[Selecta,],MN[Selecta,])


#To free some memory:
rm(Mo);rm(Di);rm(Tri);rm(Quadri);rm(Mo21);rm(Di21);rm(Tri21);rm(Quadri21);

#Safety save of input file
NormInputTrain=InputTrain



# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --
#Compute the genome wide correaltion between preditors and MNase track 
#Aim : reduce the number of predictors
# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --

A=cor(InputTrain[,1:c(dim(InputTrain)[2]-1)],InputTrain[,which(names(InputTrain)=="MN[Selecta, ]")])
A=A[order(A),]

write.table(A,file="Correlation_MNase_40_65bp_trained_set_Tacidophilum.txt")
# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --
#Plot the top and bottom 20 better correalted predictors
# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --

pdf(file="Pearson_correlation_genome_wide_top_bottom_r2d2_40_65bp.pdf",width=12,height=8)
par(mfrow=c(1,1))
barplot(A[c(1:20,(length(A)-20):length(A))],las=2,main="correlation with MNase coverage")
dev.off()

Aspear=cor(InputTrain[,1:c(dim(InputTrain)[2]-1)],InputTrain[,which(names(InputTrain)=="MN[Selecta, ]")],method="spearman")
Aspear=Aspear[order(Aspear),]

write.table(Aspear,file="Spearman_correlation_MNase_40_65bp_trained_set_Tacidophilum.txt")


# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --
#Plot the top and bottom 20 better correalted predictors
# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --

pdf(file="Spearman_correlation_genome_wide_top_bottom_r2d2_40_65bp.pdf",width=12,height=8)
par(mfrow=c(1,1))
barplot(Aspear[c(1:20,(length(Aspear)-20):length(Aspear))],las=2,main="correlation with MNase coverage")
dev.off()

# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --
#Export selected predictors ( based on correaltion absolute value) and MNase
# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --# --

#Normalize all predictors :
NormInputTrain <- NormInputTrain %>% mutate_at(-ncol(.), scale)


write.table(NormInputTrain[,which(names(NormInputTrain)%in%names(A)[c(1:40,(length(A)-39):length(A))])],file="Normalized_Input_LASSO_MOdeling_80_predictors_Two_windowSizes_30bp.txt",sep="\t",row.names=F)  


NormInputTrainExport=NormInputTrain[,which(names(NormInputTrain)%in%c(names(A)[c(1:40,(length(A)-39):length(A))],"MN[Selecta, ]"))]


names(NormInputTrainExport)[dim(NormInputTrainExport)[2]]="MNaseData"


writeMat("Normalized_MNase_and_predictors_Tacidophilum_40_65bp.mat",NormInputTrainExport=NormInputTrainExport)

