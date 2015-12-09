library(dplyr)
# Read the raw_data.tsv file
raw_data <- read.table("raw_data.tsv", header = T, fill = T)

# Make sure the data was imported in a dplyr compatible data frame
raw <- as_data_frame(raw_data)
raw %>% class() 
raw %>% sapply(class) # get the class of each column


