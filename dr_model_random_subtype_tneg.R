#### This file uses a Cox Proportional-Hazards model (no regularization or penalty) to predict patient survival (using survival analysis) using 165 drug gene signatures (i.e., RNA expression of corresponding genes in BRCA tumors)

## Output file contains three fields for each of the 165 drugs from the small molecule suite: 
# a) Drug ID - index of the drug in HMS LINCS
# b) c-index (concordance index) as a measure of model performance
# c) p-value of model performance

# install.packages("caret", dependencies = c("Depends", "Suggests"))
# source("https://bioconductor.org/biocLite.R")
# biocLite("survcomp") # for concordance.index
require(data.table) # for fread
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
library(survcomp) # for concordance.index function
library(survival) # for survival analysis
library(caret)
library(modelr) # crossv_kfold for cross validation
library(purrr) # for map, map2 functions
library(broom) # augment
# library(ROCR)
# library(pROC)
# library(limma) # quantiles normalization
# library(glmnet) # for coxnet to regularize coxph model to select relevant covariates when data is sparse


## set working directory
# setwd('/Users/savantidiot/Mangalath/Harvard/research/lsp_artem/coxanalysis/')

## Load expression data, clinical characteristics (event - death, time - days), and drug-gene signatures
load('drugrepoData.RData')

# if you want to check which data tables were loaded
# check_data <- load('drugrepoData.RData')
# check_data

## Get command line arguments
args = commandArgs(trailingOnly=TRUE)
## test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("Input file needs three arguments (start, end, step)", call.=FALSE)
} else {
  len_start = as.numeric(args[1])
  len_end = as.numeric(args[2])
  step = as.numeric(args[3])
}

# SETUP MODEL DATA: transpose data table - samples as rows, genes as columns
brca_exp <- select(brca_cancer, -one_of('gene_name'))

# remove genes where mean expression of any gene is less that 5.  Doing this avoids issues with prediction later on as we will end up with genes that cannot be tested.
brca_low_rowsum <- which(rowMeans(brca_exp) < 5) # 3801 of the 20533 genes
brca_exp <- brca_exp[-brca_low_rowsum,]
brca_cancer <- brca_cancer[-brca_low_rowsum, ] # only 41 of the genes associated with any of the drugs are part of this 3801 genes

# setup model data to contain transpose of the brca data frame with patient samples as rows and genes in columns
model_days_new <- data.table( t(brca_exp) ) # transpose only expression values, drop gene_name column
colnames(model_days_new) <- brca_cancer$gene_name  # re-assign gene names to columns (lost during transposition)

# add clinical variables to the data table
model_days_new <- model_days_new %>% mutate(days=brca_days$days, class = brca_days$status)

## Remove all patients where days count is 0 as COXNET model does not handle 0 values
# model_days_new <- model_days_new %>% filter(days != 0)

rownames(model_days_new) <- colnames(brca_exp)  # re-assign patient IDs to data table

# Separate brca data by subtype for analysis
# er_patients <- subtype_cancer %>% filter(er_status_by_ihc=="Positive") %>% select(bcr_patient_barcode)
# her2_patients <- subtype_cancer %>% filter(her2_status_by_ihc=="Positive") %>% select(bcr_patient_barcode)
tripneg_patients <- subtype_cancer %>% filter(mol_type=="TripleNegative") %>% select(bcr_patient_barcode)

## Final model data to use in training and prediction
model_data_subtype <- model_days_new[which( substr(colnames(brca_cancer),1,12) %in% tripneg_patients$bcr_patient_barcode ),]

## Setup function to get concordance index (performance metric) for random gene signatures
get_index <- function(model_matrix, num_folds){
  ## use crossv_kfold from modelr package to run n-fold cross validation
  samp_model <- crossv_kfold(model_matrix, k=num_folds, id="fold")
  
  ## train and test model on n-folds; extract average estimate of final c.index score
  # train - coxph
  samp_model <- samp_model %>% mutate(model = map(train, ~ coxph(Surv(as.numeric(days),class) ~ ., data= .)))
  # train - save training model concordance
  samp_model <- samp_model %>% mutate(train.index = map(model, ~ .x$concordance[6]))
  mean_train_cindex <- mean(unlist(samp_model$train.index))
  sd_train_cindex <- sd(unlist(samp_model$train.index))
  
  # test - predict
  samp_model <- samp_model %>% mutate(predicted = map2(model, test, ~ predict(.x, newdata=.y)))
  # get concordance index based on predicted values
  samp_model <- samp_model %>% mutate(concordance = map2(test, predicted,
                                                         ~ concordance.index(.y, surv.time = as.data.frame(.x)$days,
                                                                             surv.event = as.data.frame(.x)$class,
                                                                             method = "noether")))
  # get prediction concordance index
  samp_model <- samp_model %>% mutate(c.index = map(concordance, ~ .x$c.index))
  # get mean c.index
  mean_cindex <- mean(unlist(samp_model$c.index))
  sd_cindex <- sd(unlist(samp_model$c.index))
  indices <- list(train_index = mean_train_cindex, sd_train = sd_train_cindex,
                  test_index = mean_cindex, sd_test = sd_cindex)
  return(indices)
}

## Train model on a subset of randomly selected gene expression data and predict on testing data using the same random gene signatures
## For gene signatures up to length 200, generate 1000 randomly selected gene signatures of lengths {5, 10, 15, ....100}
## For gene signatures greater than length 100, we won't generate any random predictions for now

# generate 1000 random predictions for gene lengths of X {5, 10, 15, ....100}
# gene_len <- seq(10, 10, by = 5)
gene_len <- seq(len_start, len_end, by = step)

for (l in gene_len){
  
  ## Re-set the attributes to be captured for gene signatures of every length
  test_cindex <- vector()
  testindex_sd <- vector()
  train_cindex <- vector()
  trainindex_sd <- vector()
  drug_num <- vector()
  num_genes <- vector()
  
  for (x in 1:3){
    
    ## Set seed for reproducibility
    set.seed(x+l+713)
    
    # SETUP PREDICTORS: select sample predictor genes based on target gene set for specific drugs
    random_genes <- sample(1:ncol(model_data_subtype %>% select(-one_of("days", "class"))), size = l)
    random_matrix <- model_data_subtype %>% select(c(random_genes,"days", "class"))
    rownames(random_matrix) <- rownames(model_data_subtype)
    
    ## Print iteration number
    print(paste("Iteration number",x, "length", l))
    
    # remove gene columns where majority (>90%) of values are 0
    zero_cols <- which( colSums(select(random_matrix, -one_of("days", "class"))==0) > 0.9 * nrow(random_matrix) )
    ifelse(length(zero_cols) > 0, random_matrix <- random_matrix[,-zero_cols], random_matrix <- random_matrix)
    predictor_matrix <- random_matrix
    
    
    ## Occasionally this statement errors with the following issue: NA/NaN/Inf in foreign function call (arg 6)
    ## To prevent the error in a single iteration to break an entire run, we catch the error using try
    try_err <- try (
      model_cox <- get_index(predictor_matrix, 5) # for tneg, because of small sample size, we use 5-fold CV
      , silent = T)
    
    # skip the current iteration if there is an error, else continues
    if (class(try_err) == "try-error"){
      print("error")
      next
    }
    
    ## Store train and test cindex
    test_cindex <- c(test_cindex, model_cox$test_index)
    train_cindex <- c(train_cindex, model_cox$train_index)
    testindex_sd <- c(testindex_sd, model_cox$sd_test)
    trainindex_sd <- c(trainindex_sd, model_cox$sd_train)
    drug_num <- c(drug_num, paste0(l,x))
    num_genes <- c(num_genes, ncol(predictor_matrix)-2)
  }
  
  cox_random_subtype_df <- data.table(drug=drug_num, test.cindex=test_cindex, train.cindex=train_cindex, 
                                      testindex.sd = testindex_sd, trainindex.sd = trainindex_sd, genecount=num_genes)
  
  # write to file
  cox_random_file_subtype <- paste0("cox_random_results_subtype_tneg",l,".txt")
  write.table(cox_random_subtype_df, file = cox_random_file_subtype, row.names = F)
  
  # print(paste("File number", l, "written"))
}

