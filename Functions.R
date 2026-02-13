### ~ THE Functions ~ ###

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
library(mice)
library(missForest)
library(impute) #BioConductor
library(imputeR)
library(naniar)
library(na.tools)
library(pcaMethods)
library(missMDA)


### SET SEED
set.seed(231304)

### Defining a global set of genes that we expect to find based on results from the study of 
GenesExpected = c("Lcn12", "Nlrp12", "ENSMUSG00000084349", "ENSMUSG00000089698", "Mup8", "ENSMUSG00000081684", "ENSMUSG000000110453", "Sdk2", "Cyp2d11", "ENSMUSG00000078965", "ENSMUSG00000106576")

### DEFINE FUNCTIONS FOR DATA PROCESSING

#Remove zero rows
only_zeros_or_na <- function(row) {
  all(is.na(row) | row == 0)
}

#data preprocess function to clean and filter the metadata and countdata
pre_process <- function(data, meta) {
  
  #Remove zero rows
  data <- data %>%
    filter(!apply(., 1, only_zeros_or_na))
  
  #Name the groups based on latitude
  meta$group = "North"
  meta$group[meta$Latitude < 35] = "South"
  
  #Keep only data from extreme latitude
  meta <- rbind(meta[meta$Latitude <35,],
                meta[meta$Latitude >40,])
  
  data = data[colnames(data) %in% meta$RNAseq_Run_NCBI]
  
  return(list(data, meta))
}

#Dimensionality Reduction with Principal Component Analysis
PCA <- function(countdata, meta) {
  data <- DGEList(counts=countdata, genes=rownames(countdata), samples=meta)
  #Normlaizing the data
  dataNorm <- calcNormFactors(data)
  dataCpm <- cpm(dataNorm)
  #Get PCA Fit
  pcaFit <- pca(log2(t(dataCpm) + 0.25))
  return(scores(pcaFit))
}

#Get the Significant Genes Indicated in the Volcano Plots
get_significant_genes <- function(results) {
  new_result <- results[(results$log2FoldChange > 1 | results$log2FoldChange < 1), ]
  new_result <- new_result[log10(new_result$padj) < -5, ]
  return(new_result)
}

#Create a Nested List of Genes discovered by each algorithm
compile_gene_list <- function(resultList, nameList) {
  geneList <- list()
  for (thisResult in resultList) {
    #returns a dataframe subset of significant genes for each result from each algorithm tested
    significant_genes_df <- get_significant_genes(thisResult)
    #Create a list for the significant gene's ensemble ID's
    significant_genes <- significant_genes_df$Name
    #Append result to list to be returned
    geneList <- append(geneList, list(significant_genes))
  }
  
  #Add Originally Expected genes into the list before creating Uplot
  geneList <- append(geneList, list(GenesExpected))
  
  #Set the names
  names(geneList) <- nameList
  
  #return nested list of significant differentially expressed genes
  return(geneList)
}

### DIFFERENTIAL GENE EXPRESSION ANALYSIS

# Differential expression analysis based on the data from Mack et al. 2018
## A. Balard, 2024
## Source: https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html

#Extracting significant differential expressed genes
DEA <- function(countdata, coldata) {
  
  #Create a DESeq object
  dds <- DESeqDataSetFromMatrix(countData = round(countdata),
                                colData = coldata,
                                design = ~ Sex + group)
  ## Run analysis
  dds <- DESeq(dds)
  res <- results(dds)
  
  ## Summarize results
  summary(res)
  
  ## Convert Ensemble names to gene names (if they are known)
  results_df = as.data.frame(res)
  results_df$Ensembl_id = row.names(results_df)
  results_df = results_df[order(results_df$padj),]
  
  results_genes = gconvert(row.names(res), organism = "mmusculus",
                           target = "ENTREZGENE_ACC", filter_na = FALSE)
  # add the gene names
  results_df = merge(results_df,
                     results_genes[,c("input", "target", "name", "description")],
                     by.x = "Ensembl_id", by.y = "input")
  
  results_df$Name <- ifelse(is.na(results_df$name), results_df$Ensembl_id, results_df$name)
  
  #Return subset the results to keep only significant genes
  return(results_df[results_df$padj < 0.05 & !is.na(results_df$padj),])
  
}

### DEFINE FUNCTIONS FOR VISUALIZATIONS AND EVALUATIONS

#Volcano Plot
v_plot <- function(data) {
  EnhancedVolcano(data,
                  lab = data$Name,
                  x = 'log2FoldChange',
                  y = 'padj',
                  drawConnectors = TRUE)
  
}

# U Plot
create_UPlot <- function(geneList) {
  upset(fropcast(geneList))
}

#Function to conduct Kolmogorov-Smirnov Test and plot it 
ks_test <- function(imputed, original) {
  
  #BOTH Imputed and Original Datasets require to be transposed so that genes are treated as features (columns).
  imputed <- as.data.frame(t(imputed))
  original <- as.data.frame(t(original))
  
  # Initialize a list to store the results
  ks_test_results <- list()
  p_values <- numeric(ncol(imputed))
  column_names <- character(ncol(imputed))
  
  # Iterate over all columns
  for (col_index in seq_along(imputed)) {
    # Extract the columns, omitting missing values
    col1 <- original[[col_index]][!is.na(original[[col_index]])]
    col2 <- imputed[[col_index]]
    
    # Perform the Kolmogorov-Smirnov test
    test_result <- ks.test(col1, col2)
    
    # Store the result
    ks_test_results[[names(imputed)[col_index]]] <- test_result
    p_values[col_index] <- test_result$p.value
    column_names[col_index] <- names(imputed)[col_index]
  }
  
  # Create a data frame for plotting
  results_df <- data.frame(
    Column = column_names,
    PValue = p_values,
    Significance = ifelse(p_values < 0.05, "Significant", "Non-Significant")
  )
  
  # Calculate percentages of significant and non-significant results
  significant_count <- sum(results_df$Significance == "Significant")
  non_significant_count <- sum(results_df$Significance == "Non-Significant")
  total_count <- nrow(results_df)
  
  significant_percent <- (significant_count / total_count) * 100
  non_significant_percent <- (non_significant_count / total_count) * 100
  
  # Create a string for annotation
  annotation_text <- paste("Significant:", round(significant_percent, 2), "%, ",
                           "Non-Significant:", round(non_significant_percent, 2), "%")
  
  # Plot the p-values without a logarithmic scale and add text annotation
  ggplot(results_df, aes(x = Column, y = PValue, color = Significance)) +
    geom_point() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    geom_text(aes(x = Inf, y = Inf, label = annotation_text), 
              vjust = 2, hjust = 1.1, size = 5, color = "black") +
    labs(title = "Kolmogorov-Smirnov Test P-Values",
         x = "Column",
         y = "P-Value") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom")
}

#Variance Analysis with PCA
VariancePCAGraph <- function(med, mforest, bpca, knn, svd, pca, mean, em, cc, meta) {
  
  #Function that returns the PCA scores for each imputation technique
  PCAmed <- PCA(med, meta)
  PCAmean <- PCA(mean, meta)
  PCApca <- PCA(pca, meta)
  PCAcc <- PCA(cc, meta)
  PCAknn <- PCA(knn, meta)
  PCAbpca <- PCA(bpca, meta)
  PCAsvd <- PCA(svd, meta)
  PCAem <- PCA(em, meta)
  PCAmforest <- PCA(mforest, meta)
  
  #Plot the first two principal components of all the tested methods to summarize
  ggplot() +
    geom_point(PCAknn, mapping=aes(PC1, PC2, color="KNN")) +
    geom_point(PCAmed, mapping=aes(PC1, PC2, color="Median")) +
    geom_point(PCAmean, mapping=aes(PC1, PC2, color="Mean")) +
    geom_point(PCAbpca, mapping=aes(PC1, PC2, color="BPCA")) +
    geom_point(PCAsvd, mapping=aes(PC1, PC2, color="SVD")) +
    geom_point(PCAem, mapping=aes(PC1, PC2, color="EM")) +
    geom_point(PCApca, mapping=aes(PC1, PC2, color="PCA")) +
    geom_point(PCAmforest, mapping=aes(PC1, PC2, color="Miss Forest")) +
    geom_point(PCAcc, mapping=aes(PC1, PC2, color="Complete Cases")) +
    geom_text() 
  
}


# Plot all the recognized genes from DEA analysis

plot_barPlot <- function(data) {
  # Calculate the count of genes occurred in each sample (column-wise)
  total_gene_counts <- colSums(data)
  
  # Calculate the count of genes in each sample that are also found in Sample 1
  genes_in_original <- rownames(data[(data$Original == 1), ])
  common_gene_counts <- colSums(data[-1, , drop = FALSE])
  
  # Make sure common_gene_counts has the same length as total_gene_counts by assigning 0 for missing values
  common_gene_counts <- sapply(1:ncol(data), function(i) {
    count <- 0
    for (gene in rownames(data)) {
      if (gene %in% genes_in_original && data[gene, i] == 1) {
        count <- count + 1
      }
    }
    return(count)
  })
  
  # Create a dataframe for plotting
  plot_data <- data.frame(Method = colnames(data),
                          TotalGenes = total_gene_counts,
                          ExpectedGenes = common_gene_counts)
  
  # Melt the dataframe for stacked bar plot
  plot_data_melted <- melt(plot_data, id.vars = "Method")
  
  # Calculate position for the text labels
  plot_data_melted <- plot_data_melted %>%
    group_by(Method) %>%
    mutate(pos = cumsum(value) - 0.5 * value)
  
  #Removing Original value
  plot_data_melted <- plot_data_melted[plot_data_melted$Method != "Original", ]
  
  
  # Plot the overlapping bar plot using ggplot2
  ggplot(plot_data_melted, aes(x = Method, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.5) +
    geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5, color = "black", size = 3.5) +
    scale_fill_manual(values = c("TotalGenes" = "turquoise", "ExpectedGenes" = "plum"), 
                      labels = c("Total Genes", "Expected Genes")) +
    labs(title = "Genes found from DEA Analysis",
         x = "Imputation Methods",
         y = "Number of Genes Found",
         fill = "Legend") +
    theme_minimal()
}


#Visualize Overall DEA Results as a barplot
visualize_DEA <- function(gList) {
  
  #uncomment this line for 10% and 30% datasets
  #removing mean imputation from results
  geneList <- gList[setdiff(names(gList), "Mean")]
  
  # Convert the named nested list into a long-format dataframe
  g <- enframe(geneList, name = "list_name", value = "Genes") %>% unnest(Genes)
  
  print(g)
  # Create the binary dataframe
  df <- g %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = list_name, values_from = present, values_fill = list(present = 0))
  
  data <- as.data.frame(df)
  #use values from first column as row names
  names <- data[,1]
  data <- data[,-1]
  rownames(data) <- names
  plot_barPlot(data)
}


### DEFINE FUCNTIONS FOR IMPUTATION ALGORITHMS ATTEMPTED

# Imputation for Miss Forest Algorithm
impute_MissForest <- function(data) {
  imputed_data10 <- missForest(data, verbose = TRUE, maxiter = 10, ntree = 100, parallelize = 'variables')
  # Extract the imputed dataset
  imputed <- imputed_data10$ximp
  return(as.data.frame(t(data)))
}

# Imputation for Bayesian Principal Analysis
impute_BPCA <- function(data) {
  bpca_fit <- pca(data, method = "bpca", nPcs = 6)
  df_imputed <- completeObs(bpca_fit)
  return(as.data.frame(t(data)))
}

# Imputation with KNN with impute.knn from the Bioconductor Package
impute_knn <- function(data) {
  returned <- impute.knn(data, k=2)
  data <- returned$data
  return(as.data.frame(t(data)))
}

# Simple Statistical Median Imputation
impute_median <- function(data) {
  mImpute <- preProcess(data, method="medianImpute")
  data <- predict(mImpute, data)
  return(as.data.frame(t(data)))
}
# Simple Statistical Mean Imputation
impute_mean <- function(data) {
  for(i in 1:ncol(data)) {
    data[ , i][is.na(data[ , i])] <- mean(data[ , i], na.rm=TRUE)
  }
  return(as.data.frame(t(data)))
}

# # Singular Value Decomposition approximation estimated by the Expectation-Maximization algorithm.
impute_svd <- function(data) {
  returned <- impute_missing_val(data, algorithm = "EM-SVD", max_iter = 20000)
  data <- returned$.data
  return(as.data.frame(t(data)))
}

# The Expectation-Maximization algorithm is implemented for fitting AMMI statistical model
impute_em <- function(data) {
  returned <- impute_missing_val(data, algorithm = "EM-AMMI", max_iter = 20000)
  data <- returned$.data
  return(as.data.frame(t(data)))
}

# Function for Principal Component Analysis Imputation with package missMDA
impute_pca <- function(data) {
  imputed_data <- imputePCA(data, ncp = 5)
  data <- imputed_data$completeObs
  return(as.data.frame(t(data)))
}
