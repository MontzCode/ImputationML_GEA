### ~ THE VISUALIZATION SCRIPT ~ ###

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

### EVALUATING IMPUTATION TECHNIQUES

# LOAD DATASETS and set all negative values to 0

#Note : When Importing the Imputed Dataset, make sure you're importing the dataset that corresponds to 
#the missing value percentage you intend to analyse with the correct formatting at the end

#Loading KNN Imputed Dataset
imputedData_KNN <- read.csv("ImputedDatasets/Imputed_KNN_30%.csv", header = T, row.names = 1)
imputedData_KNN[imputedData_KNN < 0] <- 0 

#Loading SVD Imputed Dataset
imputedData_SVD <- read.csv("ImputedDatasets/Imputed_SVD_30%.csv", header = T, row.names = 1)
imputedData_SVD[imputedData_SVD < 0] <- 0 

#Loading EM-AMMI Imputed Dataset
imputedData_EM <- read.csv("ImputedDatasets/Imputed_EM_30%.csv", header = T, row.names = 1)
imputedData_EM[imputedData_EM < 0] <- 0 

#Loading Median Imputed Dataset
imputedData_MED <- read.csv("ImputedDatasets/Imputed_Median_30%.csv", header = T, row.names = 1)
imputedData_MED[imputedData_MED < 0] <- 0 

#Loading Mean Imputed Dataset
imputedData_MEAN <- read.csv("ImputedDatasets/Imputed_Mean_30%.csv", header = T, row.names = 1)
imputedData_MEAN[imputedData_MEAN < 0] <- 0 

##Loading BPCA Imputed Dataset
imputedData_BPCA <- read.csv("ImputedDatasets/Imputed_BPCA_30%.csv", header = T, row.names = 1)
imputedData_BPCA[imputedData_BPCA < 0] <- 0 

#Loading pca Imputed Dataset
imputedData_pca <- read.csv("ImputedDatasets/Imputed_pca_30%.csv", header = T, row.names = 1)
imputedData_pca[imputedData_pca < 0] <- 0 

##Loading Miss Forest Imputed Dataset
imputedData_MForest <- read.csv("ImputedDatasets/Imputed_MissForest_30%.csv", header = T, row.names = 1)
imputedData_MForest[imputedData_MForest < 0] <- 0 


# Kolmogorov-Smirnov Test to Visualize Distribution for a sample over all the genes

ks_test(imputedData_MED, countdata)
ks_test(imputedData_MEAN, countdata)
ks_test(imputedData_KNN, countdata)
ks_test(imputedData_BPCA, countdata)
ks_test(imputedData_SVD, countdata)
ks_test(imputedData_EM, countdata)
ks_test(imputedData_MForest, countdata)
ks_test(imputedData_pca, countdata)

## VARIANCE ANALYSIS

VariancePCAGraph(mforest=imputedData_MForest, med=imputedData_MED, bpca=imputedData_BPCA, knn=imputedData_KNN, svd=imputedData_SVD, pca=imputedData_pca, mean=imputedData_MEAN, em=imputedData_EM, cc=countdata_cc, meta=coldata) 

#######################################################################################################################################################################################################

#Read all the results saved 

resultsMED <- read.csv("Results/DEA_Results/results_Median_30%.csv", header = T)

resultsMEAN <- read.csv("Results/DEA_Results/results_Mean_30%.csv", header = T)

resultsKNN <- read.csv("Results/DEA_Results/results_KNN_30%.csv", header = T)

resultsSVD <- read.csv("Results/DEA_Results/results_SVD_30%.csv", header = T)

resultsEM <- read.csv("Results/DEA_Results/results_EM_30%.csv", header = T)

resultsMForest <- read.csv("Results/DEA_Results/results_MissForest_30%.csv", header = T)

resultsBPCA <- read.csv("Results/DEA_Results/results_BPCA_30%.csv", header = T)

resultspca <- read.csv("Results/DEA_Results/results_pca_30%.csv", header = T)


# Club all significant results from DEA into a single list for Visualization Graphs
results <- list(resultsMED, resultsMEAN, resultsKNN, resultsSVD, resultsEM, resultsMForest, resultsBPCA, resultspca)

# We need to save the names in order of results df
nameList <- list("Median", "Mean", "KNN", "SVD", "EM", "Miss_Forest", "BPCA", "pca", "Original")

# Compile a Gene List
geneList <- compile_gene_list(results, nameList)

#Visualize DEA with VOLCANO for all the results individually

v_plot(resultsMED)
v_plot(resultsMEAN)
v_plot(resultsKNN)
v_plot(resultsSVD)
v_plot(resultsEM)
v_plot(resultsMForest)
v_plot(resultspca)
v_plot(resultsBPCA)

#Visualize Overall DEA Results as a barplot
visualize_DEA(geneList)
