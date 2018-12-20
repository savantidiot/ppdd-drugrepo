#### This file executes the Cox-ph regression algorithm on 165 drug gene signatures and stores the output in a text file.
## Output file contains three fields: a) Gene ID, b) c-index as a measure of model performance, c) p-value of model performance
## Note: Please download the source TCGA data from the following link prior to running the code:  https://www.synapse.org/#!Synapse:syn4976369
## If you are unable to download TCGA data, please contact author at pradeep underscore mangalath at hms dot harvard dot edu
## Source files needed include: a) TCGA data (clinical and expression), b) Drug gene signatures from .gmt file available on the Git repo


# install.packages("caret", dependencies = c("Depends", "Suggests"))
# source("https://bioconductor.org/biocLite.R")
# biocLite("survcomp") # for concordance.index
library(survcomp)
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


# filter patients total survival in days (for patients in class 'Alive', take value from "days_to_last_followup"; for patients in class 'Alive', take value from "days_to_death")
brca_days <- brca_clin_data %>% filter( bcr_patient_barcode %in% substr(colnames(brca_data_cancer), 1, 12) ) %>% select(c(bcr_patient_barcode,as.numeric(days_to_last_followup), as.numeric(days_to_death)))
# brca_days$days_to_death <- brca_days$days_to_death %>% replace_na(0)
brca_days <- brca_days %>% mutate(days = 0)
# brca_days$days[which(brca_vital=="Dead")] <- brca_days$days_to_death[which(brca_vital=="Dead")]
# brca_days$days[which(brca_vital=="Alive")] <- brca_days$days_to_death[which(brca_vital=="Alive")]
 
# brca_days <- brca_days %>% mutate(days[which(brca_vital=="Dead")] = days_to_death[which(brca_vital=="Dead")], days[which(brca_vital=="Alive")] = days_to_last_followup[which(brca_vital=="Alive")])
brca_days[which(brca_vital=="Dead"), "days"] <- brca_days[which(brca_vital=="Dead"), "days_to_death"]
brca_days[which(brca_vital=="Alive"), "days"] <- brca_days[which(brca_vital=="Alive"), "days_to_last_followup"]

# final checks on data - remove negative days (make zero) and change NA values to zero (similar to lost to follow up)
days_na <- which(is.na(brca_days$days))
ifelse(length(days_na) > 0, brca_days[days_na, "days"] <- brca_days[days_na,"days_to_last_followup"], brca_days <- brca_days)
days_lessthanzero <- which(brca_days$days < 0)
ifelse(length(days_lessthanzero) > 0, brca_days[days_lessthanzero, "days"] <- 0, brca_days <- brca_days)
brca_days$days <- as.numeric(brca_days$days)
brca_days <- brca_days %>% mutate(status = brca_vital)


# cox model requires censoring status to be 0 or 1; here we censor if Alive (0) as we're modeling surivival outcomes
brca_days$status <- ifelse(brca_days$status =="Alive",0,1)

# SETUP MODEL DATA: transpose data table - samples as rows, genes as columns
model_days <- t(brca_data_cancer)
colnames(model_days) <- tcga_data$gene_id  # re-assign gene names to columns (lost during transposition)
model_days_new <- data.table(model_days)
model_days_new <- model_days_new  %>% mutate(days=brca_days$days, class = brca_days$status)
rownames(model_days_new) <- rownames(model_days)


####### Visualizing the survival data
brca_surv <- Surv(as.numeric(model_days_new$days), model_days_new$class)


#### For 165 drug gene signatures
Cox_cindex <- vector()
Cox_pval <- vector()
set.seed(209)

for (x in 1:165){
  ## Iterate over the 165 genes
  print(paste("Iteration number",x))
  gene_list <- gene_array[[x]]
  
  # SETUP PREDICTORS: select sample predictor genes based on target gene set for specific drugs
  # random_genes <- sample(1:ncol(model_days), size = 50)
  random_genes <- grep(paste0("\\b",gene_list,"\\b",collapse="|"), colnames(model_days), ignore.case = T)
  random_matrix <- model_days_new %>% select(c(random_genes, "days", "class"))
  rownames(random_matrix) <- rownames(model_days)
  
  # remove gene columns where majority (>70%) of values are 0
  zero_cols <- which( colSums(random_matrix==0) > 0.9 * nrow(random_matrix) )
  ifelse(length(zero_cols) > 0, random_matrix <- random_matrix[,-zero_cols], random_matrix <- random_matrix)
  predictor_matrix <- random_matrix
  # predictor_matrix <- predictor_matrix %>% mutate(status=brca_vital)
  # predictor_matrix <- predictor_matrix %>% mutate(newclass=rep(1,1097))
  
  
  ### SPLIT TRAINING AND TESTING DATA: create training and test data sets
  # set.seed(103)
  modTrain <- createDataPartition(
    y = predictor_matrix$class,
    p= .7,
    list = FALSE
  )
  
  training <- predictor_matrix[modTrain,]
  testing <- predictor_matrix[-modTrain,]
  

  ### TRAIN MODEL: Using COX model
  Model_Cox <- coxph(Surv(as.numeric(days),class) ~ ., data=training)
  
  
  ### PREDICT MODEL: using predict and ROC
  Model_Cox_Pred <- predict(Model_Cox, newdata = testing)
  
  # Determine concordance
  Model_Cox_cindex = concordance.index(Model_Cox_Pred, surv.time = testing$days,
                                        surv.event=testing$class, method = "noether")
  
  ## Store AUC and ROC values in a list
  Cox_cindex <- c(Cox_cindex, Model_Cox_cindex$c.index)
  Cox_pval <- c(Cox_pval, Model_Cox_cindex$p.value)
  
}

cox_drug_random_df <- data.table(drug=gene_link,cindex=Cox_cindex, pval=Cox_pval)

# write to file
cox_drug_file <- "cox_drug_results.txt"
write.table(cox_drug_random_df, file = cox_drug_file, row.names = F)


############## End of Code #############
# ####### Visualizing the survival data
# 
# brca_surv_training <- Surv(as.numeric(training$days), training$class)
# brca_surv_training <- survfit(brca_surv_training~1, data = training)
# autoplot(brca_surv_training, conf.int = T, color="red") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
#   xlab("Time (days)") + ylab("Proportion of patients alive (percentage)") +
#   ggtitle(label = "Survival curve for 1097 patient samples based on TCGA BRCA data")
# 
# brca_surv_testing <- Surv(as.numeric(testing$days), testing$newclass)
# brca_surv_testing <- survfit(brca_surv_testing~status, data = testing)
# autoplot(brca_surv_testing, conf.int = F) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
#   xlab("Time (days)") + ylab("Proportion of patients without disease (DFS)") +
#   ggtitle(label = "Disease Free Survival Analysis of AML and ALL patients by group")
# 
# 
# clin_days <- ifelse(is.na(brca_clin_data$days_to_death), as.numeric(brca_clin_data$days_to_last_followup), as.numeric(brca_clin_data$days_to_death))


