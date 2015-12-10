#!/bin/env Rscript
#==================================================================
# Author: Yassine Souilmi
# 
# This is a response to the exercise initially posted here:
#        https://github.com/mpg-age-bioinformatics/exercise
#==================================================================

#########################
# Reading the inputs    #
#########################
#==================================================================
### Cleanup the environment before the execution
rm(list = ls())

### Read the input arguments
args<-commandArgs(TRUE)

### Make sure that the user provided an argument
if (length(args) == 0) { 
  print("No input argument was provided!")
} 

### Read the provided data-path
input=toString(args[1])
# input = '~/bin/mpg_ex/mpg/raw_data.tsv'

### Change the working-directory to the specified path
x = unlist(strsplit(input, split='/', fixed=TRUE))
wd = paste(x[1:(length(x)-1)], collapse = '/')
setwd(wd)
#==================================================================
#########################
# Loading the libraries #
#########################

### Load needed libraries 
library(dplyr)  # for data manipulation
library(ggvis)  # for visualizaiton 
#==================================================================
#######################################################
# Reading the data  and taking a sneak peek at it     #
#######################################################

# Read the raw_data.tsv file
raw_data <- read.table("x", header = T, fill = T)

# Make sure the data was imported in a dplyr compatible data frame
raw <- raw_data %>% as_data_frame()

# remove raw_data dataframe
rm(raw_data)

# Check the content of the raw dataframe
raw %>% class() 
raw %>% sapply(class) # get the class of each column
raw %>% head()

#==================================================================
#########################
#     Analysis          #
#########################

# creat counts subsets per month and isoform
august_t <- raw %>% select(August_1t:August_4t)       # August isoform_normal
december_t <- raw %>% select(December_1t:December_4t)   # December isoform_normal
august_m <- raw %>% select(August_1m:August_4m)       # August isoform_malignant
december_m <- raw %>% select(December_1m:December_4m)   # December isoform_malignant

# Visualize the data
raw %>% ggvis(~August_1t, ~August_1m)


