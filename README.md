#Antoine Hocher, Friday 1 Novembre 2019, London, UK | MRC, London School of Medical Science 
#####################################################
#Scripts associated to Hocher et al. 2019, The DNA-binding protein HTa from Thermoplasma acidophilum is an archaeal histone analog


#in order to run LASSO model prediction :

First run LASSO_Input_file_generation.R, the ouput of which should be used as input of AH_LASSO_script.m, which runs in matlab R2018a
The output of AH_LASSO_script.m (the coefficient file), should at last be used as input of LASSO_output_file_generation.R which
will produce a bigwig tracks from the LASSO coefficients and the computed Kmers abundancy


#In order to detect and score peaks on two different mnase tracks , and compute their relative asymmetry :

run Peak_detection_and_scoring_on_indep_bwFile.R, which relies on a modified version of the published Bioconductor package NucleR
(Flores O, Orozco M (2011). “nucleR: a package for non-parametric nucleosome positioning.” Bioinformatics, 27, 2149–2150. doi: 1093/bioinformatics/btr345)

#In order to use the published plot2DO function (https://github.com/rchereji/plot2DO) (freely available MIT licence) on single end reads I had to implement small modifications (all indicated in the core code)
