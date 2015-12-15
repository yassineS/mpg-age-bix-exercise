#!/bin/env Rscript

################################################################
# Author: Yassine Souilmi                                      #
#                                                              #
# This is a response to the exercise initially posted here:    #
#        https://github.com/mpg-age-bioinformatics/exercise    #
################################################################

#########################
# Reading the inputs    #
#########################

# Cleanup the environment before the execution    
rm(list = ls())

# Read the input arguments
args = commandArgs(TRUE)

# Make sure that the user provided an argument
if (length(args) == 0) { 
  print("No input argument was provided!")
  print("./scd.R [input data]")
} 

# Read the provided data-path
input=toString(args[1])
# debug: input   = '~/bin/mpg_ex/mpg/raw_data.tsv'


# Change the working-directory to the specified path
x = unlist(strsplit(input, split='/', fixed=TRUE))
wd = paste(x[1:(length(x)-1)], collapse = '/')
setwd(wd)

#########################
# Loading the libraries #
#########################

# Load needed libraries:
# dplyr: for better handle on data manipulation
# edgeR: for binomial tests
# DESeq: for building the object that holds the dataset \
# to be analysed by edgeR and the subsequent calculations \
# performed on the dataset
library(dplyr)
source("https://bioconductor.org/biocLite.R") # make sure bioconductor is available
biocLite() # install core packages
biocLite("DESeq", "edgeR") # install libraries
library("DESeq", "edgeR") # load libraries

####################
# Reading the data #
####################

# Read the raw_data.tsv file
raw_data = read.delim(input, header = T, sep = "\t", 
                       dec = ".", fill = T, 
                       na.strings = c("inf")) # get rid of the "inf" values

# Make sure the data was imported in a dplyr compatible data frame
counts = raw_data %>% as_data_frame()

# Selecting the malignant isoform data
malignant = counts %>% select(August_1m:December_4m)

# remove raw_data dataframe
rm(raw_data, args, counts, x)

#####################
# Data exploration  #
#####################
# Visualize graphically the data
png("data_viz_before.png", width = 2400, height = 500, units = "px")
counts %>% boxplot(~August_1m) # this is a bug that we can take advantage of
dev.off()

#########################
#     Data cleaning     #
#########################
# Since the read counts are non-negative integers (Anders, 2010) \
# negative values. missing values (NAs) are replaced with '0's \
# (incompatible with DGEList and edgeR)
malignant[malignant<=0]     = 0
malignant[is.na(malignant)] = 0

# Visualize graphically the data to make sure the cleaning worked
png("data_viz_after.png", width = 2400, height = 500, units = "px")
malignant %>% boxplot(~August_1m) # this is a bug that we can take advantage of
dev.off()

####################
#     Analysis     #
####################
# 1. Make a DGEList for edgeR:
# 1.1. Define the groups
malignant_counts_groups = c(rep("August_m",4),rep("December_m",4))

# 1.2. Create the DGEList object
malignant_DGEList = malignant %>% 
  DGEList(group=malignant_counts_groups)

# 2. Run edgeR
malignant_DGEList = malignant_DGEList %>% 
  calcNormFactors() %>% 
  estimateCommonDisp(verbose=T) %>% 
  estimateTagwiseDisp() 

malignant_tgw = malignant_DGEList %>% exactTest()

# 3. Summarize the findings
# Classify a series of related differential expression statistics \
# as up, down or not significant
dt = decideTestsDGE(malignant_tgw, p.value=0.01, adjust.method = "fdr")


# Get the 'fdr' adjusting the pvalue \
# using Benjamini & Hochberg (1995) ("BH" or its alias "fdr") method
PValue_fdr <- p.adjust(method="fdr",p=malignant_tgw$table$PValue)

# Create a new dataframe to help summarize the findings
ind = c(paste("ind", 1:dim(malignant)[1], sep = '-'))
results = data.frame(ind, malignant_tgw$table, 
                     "twd"=malignant_DGEList$tagwise.dispersion,
                     "PValue_fdr"=PValue_fdr,
                     "Decide_test"=dt,
                     malignant)

results_summary = results %>% filter(PValue_fdr<=0.01)

#######################
# Saving the results  #
#######################
# Write out the results to a dataframe
write.table(results, file="results.tsv", quote=F)
write.table(results_summary, file="results_summary.tsv", quote=F)
sessionInfo()
