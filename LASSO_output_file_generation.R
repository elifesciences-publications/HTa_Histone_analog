#Antoine Hocher MRC, 2019
rm(list=ls())

##########################################################
#Compute bigwig file from Coefficient computed in matlab : 
##########################################################

ComputeBigWig <- function(CoefsFile,MNFile,ChromosomeName,Note){
  Coefs=read.table(CoefsFile,header=T,sep="\t")
  
  Coefs=Coefs[which(Coefs$Value!=0),]
  WholeG=cbind(Mo,Di,Tri,Quadri,Mo21,Di21,Tri21,Quadri21)
  
  
  WholeG=WholeG[,which(names(WholeG)%in%as.character(Coefs$Name[Coefs$Name!="Intercept"]))]
  WholeG[,1:dim(WholeG)[2]]=scale(WholeG[,1:dim(WholeG)[2]])
  
  
  Coef_noInt=Coefs[which(Coefs$Name != "Intercept"),]
  
  WholeG$PredictMN=Coefs[which(Coefs$Name == "Intercept"),]$Value
  
  #Computing the modeled data genome wide :
  for(i in 1:length(Coef_noInt$Name)){
    WholeG$PredictMN=WholeG$PredictMN+Coef_noInt[i,]$Value*WholeG[,which(names(WholeG)==Coef_noInt[i,]$Name)]
  }
  
  #Creating a table out of this data :
  Predicted_val_Tab=as.data.frame(matrix(nrow=dim(WholeG)[1],ncol=0))
  Predicted_val_Tab$Chromosome="ChromosomeName"
  Predicted_val_Tab$Start=1:dim(WholeG)[1]
  Predicted_val_Tab$End=1:dim(WholeG)[1]
  Predicted_val_Tab$score=WholeG$PredictMN
  BABO=makeGRangesFromDataFrame(Predicted_val_Tab)
  score(BABO)=WholeG$PredictMN
  
  seqlengths(BABO)=dim(WholeG)[1]
  
  
  export.bw(BABO,paste(gsub(pattern = "\\.bw$", "", basename(MNFile)),"_Predicted_",Note,"80_pred.bw",sep=""))
}



ComputeBigWig("/Users/ahocher/Documents/TEAC_HTa_Study/Results/Modeling_HTa_binding/Fitting_Coefficients_Output/Coefficients_80_pred_LASSO_r2d2_t15_40_65bp.txt","/Users/ahocher/Documents/TEAC_HTa_Study/Results/MNase_Experiments/Sequencing/2018_07_PE_TEAC/Coverage_tracks_Merged_reads_Bowtie2/RGC_Input_Zscore_Normalized/40_65/Merged_Trimmed_r2_d2_t15_GTGAAA_L001_R_001_srt_RPGC_45_60bpNormalizedNormalized_by_Input_top_decileInput_normalized_Zscore.bw","")