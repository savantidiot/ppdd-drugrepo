# ppdd-drugrepo
Code for drug repurposing methods used in PPDD course

This section contains code and supporting files used to generate the analysis and results for the Principles and Practices of Drug Development course at MIT.

Code includes:

- pre-processing of BRCA data
- pre-processing of HMS-LINCS data
- implementation of Random Forests algorithm to classify tumor status (sample_analysis_RF.r)
- implementation of Logistic Regression model to classify survival status (sample_analysis_LOGIT.r)
- implementation of Cox-PH Regression algorithm to predict survival time (sample_analysis_COX.R)


Data includes: 
- TCGA tissue sample data available to download at https://www.synapse.org/#!Synapse:syn4976369
- HMS-LINCS drug-gene signature data (*.gmt file)
