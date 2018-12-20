#### This file executes the Logistic Regression algorithm on 165 drug gene signatures and stores the output in a text file.
## Output file contains five fields: a) Gene ID, b) ROC from training model ,c) Predicted AUC score, d) Accuracy of prediction, e) p-value of NIR
## Note: Please download the source TCGA data from the following link prior to running the code:  https://www.synapse.org/#!Synapse:syn4976369
## If you are unable to download TCGA data, please contact author at pradeep underscore mangalath at hms dot harvard dot edu
## Source files needed include: a) TCGA data (clinical and expression), b) Drug gene signatures from .gmt file available on the Git repo

# install.packages("caret", dependencies = c("Depends", "Suggests"))
require(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(edgeR) # for differential gene analysis
library(survival) # for survival analysis
library(broom) # for tidy()
library(tidyr)
library(caret)
library(ROCR)
library(pROC)


## Parses a .gmt file and puts it into the list format
n_file <- 'Nienke_10genes.gmt'
# Get gene list for drug from Nienke's file
n_data <- fread(n_file, header=F, sep='\t', check.names = FALSE, nrows=165, fill = TRUE)
gene_array <- list()
gene_link <- vector()
for (r in 1:nrow(n_data)){
  gene_temp <- as.matrix(n_data[r,3:ncol(n_data),with=F])
  gene_temp <- gene_temp[ which(gene_temp != "") ]
  g_link <- strsplit(as.character(n_data[r,1]),"=")[[1]][2] # split the first term, i.e., link ("http://lincs.hms.harvard.edu/db/?search=10001") to obtain ID
  gene_array[r] <- list(gene_temp)
  gene_link <- c(gene_link, g_link)
}


# Get TCGA and clinical data from raw files
# setwd('/Users/savantidiot/Mangalath/Harvard/research/lsp_artem/tcga_data')
clin_file = 'clinical_PANCAN_patient_with_followup.tsv'
tcga_file = 'EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'
clin_data <- fread(clin_file, header=T, sep='\t', check.names = FALSE, na.strings = c("NA", "[Not Available]","[Not Applicable]"))
tcga_data <- fread(tcga_file, header=T, sep='\t', check.names = FALSE)


# get subset of brca data using values from clinical data barcodes
# tcga_new <- tcga_data  # store tcga data into new data table
brca_clin_data <- clin_data %>% filter(acronym=="BRCA")  # select only BRCA patients from clinical data set
brca_cols_tcga <- which( substr(colnames(tcga_data), 1, 12) %in%  brca_clin_data$bcr_patient_barcode )  # get index of columns in TCGA data that are BRCA related
brca_data <- tcga_data %>% select(brca_cols_tcga)  # select columns based on index for BRCA samples
rownames(brca_data) <- tcga_data$gene_id  # assign gene names to brca data


# identify tumor samples from normal samples for classification and store as factor
norm_codes <- c("11","12","13","14","15","16","17","18","19")  # In the TCGA barcode, characters 14-15 indicate whether tissue is normal or cancerous.  Characters "11"-"19" indicate normal samples
brca_normal_cols <- which( substr(colnames(brca_data), 14, 15) %in% norm_codes )  # capture index of normal samples
brca_type <- rep("cancer", ncol(brca_data))
brca_type[brca_normal_cols] <- "normal"
brca_type <- factor(brca_type)
brca_data_cancer <- brca_data %>% select(-brca_normal_cols)
brca_data_normal <- brca_data %>% select(brca_normal_cols)

# tissue samples contain a few duplicates from the same patient -- we will remove these duplicates to retain one sample per patient
brca_tiss_dup <- names(which(table( substr(colnames(brca_data_cancer), 1, 12) ) > 1))  # get counts of substring; search for counts >1; get the column names of those indexes
brca_tiss_dup <- grep(paste0("\\b",brca_tiss_dup,"\\b",collapse="|"), colnames(brca_data_cancer), ignore.case = T)  # gets all duplicates
brca_tiss_dup <- brca_tiss_dup[c(TRUE,FALSE)] # removes alternate indexes as the duplicates are adjacent to each other
brca_data_cancer <- brca_data_cancer %>% select(-brca_tiss_dup)  # remove duplicates from tissue sample data table


# filter patients vital status ("Alive" vs. "Dead") for the de-duplicated sample data
brca_vital <- brca_clin_data %>% filter( bcr_patient_barcode %in% substr(colnames(brca_data_cancer), 1, 12) ) %>% select(c(bcr_patient_barcode,vital_status))
brca_vital <- factor(brca_vital$vital_status)  # make vital status a factor type vector

# SETUP MODEL DATA: transpose data table - samples as rows, genes as columns
model_vitals <- t(brca_data_cancer)
colnames(model_vitals) <- tcga_data$gene_id  # re-assign gene names to columns (lost during transposition)
model_vitals_new <- data.table(model_vitals)
model_vitals_new <- model_vitals_new  %>% mutate(class = brca_vital)
rownames(model_vitals_new) <- rownames(model_vitals)


#### For 165 drug gene signatures
RF_Vital_modelroc <- vector()
RF_Vital_auc <- vector()
RF_Vital_accuracy <- vector()
RF_Vital_nir_pval <- vector()
set.seed(766)
for (x in 1:165){
  ## Iterate over the 165 genes
  print(paste("Iteration number",x))
  gene_list <- gene_array[[x]]
  
  # SETUP PREDICTORS: select sample predictor genes based on target gene set for specific drugs
  # random_genes <- sample(1:ncol(model_vitals), size = 50)
  random_genes <- grep(paste0("\\b",gene_list,"\\b",collapse="|"), colnames(model_vitals), ignore.case = T)
  random_matrix <- model_vitals_new %>% select(c(random_genes, "class"))
  rownames(random_matrix) <- rownames(model_vitals)
  # remove gene columns where majority (>70%) of values are 0
  zero_cols <- which( colSums(random_matrix==0) > 0.9 * nrow(random_matrix) )
  ifelse(length(zero_cols) > 0, random_matrix <- random_matrix[,-zero_cols], random_matrix <- random_matrix)
  predictor_matrix <- random_matrix
  
  
  ### SPLIT TRAINING AND TESTING DATA: create training and test data sets
  # set.seed(103)
  modTrain <- createDataPartition(
    y = predictor_matrix$class,
    p= .7,
    list = FALSE
  )
  
  training <- predictor_matrix[modTrain,]
  testing <- predictor_matrix[-modTrain,]
  
  
  ### Using SMOTE method - hybrid approach
  ctrl <- trainControl(
    method = "repeatedcv", # K fold cross validation with default of 10 fold; use number to change fold count
    repeats = 3,
    classProbs = TRUE,  # to compute predicted class probabilities
    summaryFunction = twoClassSummary, # for ROC performance
    sampling = "smote"
  )
  
  # TUNE: tune grid parameter
  tgrd <- expand.grid(.cost = 1,
                      .loss = "L2_primal",
                      .epsilon = seq(0.001, 0.01, length.out = 5))
  
  ### TRAIN MODEL: Using Random Forest
  Model_Vital_RF <- train(
    class ~ .,
    data = training,
    method = "regLogistic",
    preProc = c("center", "scale"),
    tuneGrid = tgrd,
    trControl = ctrl,
    metric = "ROC"
  )
  
  
  ### PREDICT MODEL: using predict and ROC
  Model_Vital_Pred <- predict(Model_Vital_RF, newdata = testing)
  Model_Vital_Confusion <- confusionMatrix(Model_Vital_Pred, testing$class)
  Model_Vital_Pred.prob <- predict(Model_Vital_RF, newdata = testing, type = "prob")
  Model_Vital_roc <- roc(testing$class, Model_Vital_Pred.prob$Alive)
  # plot(Model_Vital_roc, print.thres="best")
  
  ## Store AUC and ROC values in a list
  RF_Vital_modelroc <- c(RF_Vital_modelroc, Model_Vital_RF$results$ROC[1])
  RF_Vital_auc <- c(RF_Vital_auc, Model_Vital_roc$auc[1])
  RF_Vital_accuracy <- c(RF_Vital_accuracy, Model_Vital_Confusion$overall["Accuracy"])
  RF_Vital_nir_pval <- c(RF_Vital_nir_pval, Model_Vital_Confusion$overall["AccuracyPValue"])
  
}

rfvital_drug_df <- data.table(drug=gene_link,modelroc=RF_Vital_modelroc, auc=RF_Vital_auc, accuracy=RF_Vital_accuracy, nirpval=RF_Vital_nir_pval)

# write to file
rfvital_drug_file <- "vital_drug_results_actual.txt"
write.table(rfvital_drug_df, file = rfvital_drug_file, row.names = F)

