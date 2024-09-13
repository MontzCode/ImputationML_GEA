### ~ THE FINAL SCRIPT ~ ###

#packages

library(data.table)
library(dplyr)
library(tidyverse)
library(tibble)
library(magrittr)

#Parallel processing libraries
library(parallel)
library(doParallel)

#Libraries for the genomic analysis
library(BiocManager)
library(gprofiler2)
library(DESeq2)
library(caret)

#Visualization Libraries
library(ggplot2)
library(gridExtra)
library(EnhancedVolcano)
library(edgeR)
library(UpSetR)

#imputation libraries
library(metan)
library(VIM)
library(missForest)
library(mice)
library(impute) #BioConductor
library(imputeR)
library(naniar)
library(na.tools)
library(pcaMethods)
library(missMDA)

source("Functions.R")

#Import Selected Dataset and Metadata from text files
#make sure you load the dataset corresponding to the missingness degree you wish to analyze
countdata <- read.table("data/countdata30.txt", header = T, row.names = 1)
metadata <- read.table("data/metadata.txt", header = T)

#Pre-processing the Data
returned <- pre_process(countdata, metadata)
countdata <- returned[[1]] # countdata: a table with the fragment counts
coldata <- returned[[2]] # coldata: a table with information about the samples

#making sure the columns of the count matrix correspond to the rows of the sample information table.
colnames(countdata) == coldata$RNAseq_Run_NCBI
coldata <- coldata[match(colnames(countdata), coldata$RNAseq_Run_NCBI), ]
colnames(countdata) == coldata$RNAseq_Run_NCBI

##complete cases
countdata_cc = na.omit(countdata)

#UNDERSTANDING MISSINGNESS PATTERN

#Stats
#summary(countdata) #summary statistics
#sum(is.na(countdata)) #sum of all missing value in entire dataset
#colSums(is.na(countdata)) #Count missing values per sample

#MCAR test 
mcar_test_result <- mcar_test(countdata)

#Visualisation of data missingness pattern
vis_miss(countdata, warn_large_data = FALSE)

### PARALLEL COMPUTATION FOR EFFICIENCEY

# Detect the number of physical cores
physical_cores <- detectCores(logical = FALSE)
# Compute one less than the number of physical cores
num_cores <- physical_cores - 1 # use one less than the number of available cores
cl <- makeCluster(num_cores) #Create Clusters
registerDoParallel(cl) #Register Parallel

### IMPUTATION TECHNIQUES TESTED

#Note : in all fucntions the data passed is transposed and then imputed since we need the imputations to consider the genes as features
countdata <- t(countdata)

# Note : When Exporting the Imputed Dataset, make sure you're exporting the dataset that corresponds to 
#the missing value percentage you intend to analyse with the correct formatting at the end

#1 Perform Imputations with Median Imputation
imputedData_MED <- impute_median(countdata)
#write.csv(imputedData_MED, file = "ImputedDatasets/Imputed_Median_30%.csv", row.names = T, col.names = T)

#2 Perform Imputation with Singular Value Decomposition Algorithm
imputedData_SVD <- impute_svd(countdata)
#write.csv(imputedData_SVD, file = "ImputedDatasets/Imputed_SVD_30%.csv", row.names = T, col.names = T)

#3 Perform Imputation with EM Algorithm fitted on the AMMI statistical model.
imputedData_EM <- impute_em(countdata)
#write.csv(imputedData_EM, file = "ImputedDatasets/Imputed_EM_30%.csv", row.names = T, col.names = T)

#4 Perform Imputation with Clustering Based Method : K-Nearest-Neighbors
imputedData_KNN <- round(impute_knn(countdata))
#write.csv(imputedData_KNN, file = "ImputedDatasets/Imputed_KNN_30%.csv", row.names = T, col.names = T)

#5 Perform Imputation with Tree Based Method : Miss Forest
imputedData_MForest <- impute_MissForest(countdata)
#write.csv(imputedData_MForest, file = "ImputedDatasets/Imputed_MissForest_30%.csv", row.names = T, col.names = T)

#6 Perform Imputation with BPCA : combined an EM approach for PCA with a Bayesian model
imputedData_BPCA <- impute_BPCA(countdata)
#write.csv(imputedData_BPCA, file = "ImputedDatasets/Imputed_BPCA_30%.csv", row.names = T, col.names = T)

#7 perform Imputation with Maximum Likelihood Imputation from MissMDA Package
imputedData_pca <- impute_pca(countdata)
#write.csv(imputedData_pca, file = "ImputedDatasets/Imputed_pca_30%.csv", row.names = T, col.names = T)

#8 Perform Imputations with Mean Imputation
imputedData_MEAN <- impute_mean(countdata)
#write.csv(imputedData_MEAN, file = "ImputedDatasets/Imputed_Mean_30%.csv", row.names = T, col.names = T)

#Change countdata back to dataframe format
countdata <- data.frame(t(countdata))

# Shut down the cluster
stopCluster(cl)

### DIFFERENTIAL GENE EXPRESSION ANALYSIS ON IMPUTED DATASETS

# Note : When Exporting the Results Dataset, make sure you're exporting the dataset that corresponds to 
#the missing value percentage you intend to analyse with the correct formatting at the end

resultsMED <- DEA(imputedData_MED, coldata)
#write.csv(resultsMED, file = "Results/DEA_Results/results_Median_30%.csv", row.names = F)

resultsMEAN <- DEA(imputedData_MEAN, coldata)
#write.csv(resultsMEAN, file = "Results/DEA_Results/results_Mean_30%.csv", row.names = F)

resultsKNN <- DEA(imputedData_KNN, coldata)
#write.csv(resultsKNN, file = "Results/DEA_Results/results_KNN_30%.csv", row.names = F)

resultsSVD <- DEA(imputedData_SVD, coldata)
#write.csv(resultsSVD, file = "Results/DEA_Results/results_SVD_30%.csv", row.names = F)

resultsEM <- DEA(imputedData_EM, coldata)
#write.csv(resultsEM, file = "Results/DEA_Results/results_EM_30%.csv", row.names = F)

resultsMForest <- DEA(imputedData_MForest, coldata)
#write.csv(resultsMForest, file = "Results/DEA_Results/results_MissForest_30%.csv", row.names = F)

resultsBPCA <- DEA(imputedData_BPCA, coldata)
#write.csv(resultsBPCA, file = "Results/DEA_Results/results_BPCA_30%.csv", row.names = F)

resultspca <- DEA(imputedData_pca, coldata)
#write.csv(resultspca, file = "Results/DEA_Results/results_pca_30%.csv", row.names = F)
