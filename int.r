# Make sure the data was imported in a dplyr compatible data frame
raw <- raw_data %>% as_data_frame()

# remove raw_data dataframe
rm(raw_data, args)

# Check the content of the raw dataframe
raw %>% class()
raw %>% sapply(class) # get the class of each column
raw %>% head()

#################################
#     Visual screening          #

# Graphical discovery
# Visualize graphically the data
png("data_viz.png", width = 2400, height = 500, units = "px")
raw %>% boxplot(~August_1m)
dev.off()

##########################
#     Data cleaning     #

## Since the read counts are non-negative integers (Anders, 2010) \
## negative values are eleminated, and replaced with 'NA's
raw[raw<=0] <- 0
raw[is.na(raw)] <- 0

# cutoff above min count (34)

pval <- matrix(0, nrow = 16, ncol = 16)
colnames(pval) <- colnames(raw)
rownames(pval) <- colnames(raw)

#####
# Whithin: unpaired t-test
# Malignant
pval[1,2] <- t.test(raw$August_1m, raw$August_2m)$p.value
pval[1,3] <- t.test(raw$August_1m, raw$August_3m)$p.value
pval[1,4] <- t.test(raw$August_1m, raw$August_4m)$p.value
pval[2,3] <- t.test(raw$August_2m, raw$August_3m)$p.value
pval[2,4] <- t.test(raw$August_2m, raw$August_4m)$p.value
pval[3,4] <- t.test(raw$August_3m, raw$August_4m)$p.value

pval[5,6] <- t.test(raw$December_1m, raw$December_2m)$p.value
pval[5,7] <- t.test(raw$December_1m, raw$December_3m)$p.value
pval[5,8] <- t.test(raw$December_1m, raw$December_4m)$p.value
pval[6,7] <- t.test(raw$December_2m, raw$December_3m)$p.value
pval[6,8] <- t.test(raw$December_2m, raw$December_4m)$p.value
pval[7,8] <- t.test(raw$December_3m, raw$December_4m)$p.value

which(pval>0.05, arr.in=TRUE)

#####
# Selecting the malignant isoform data

malignant <- raw[, 1:8]

summary(rowSums(malignant))

# Define the groups
group <- c(rep("August_m",4),rep("December_m",4))

#create the DGEList object
d <- DGEList(counts = malignant, group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d,verbose=T)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d)
summary(decideTestsDGE(de.tgw, p.value=0.01))
topTags(de.tgw, n=10)

# Create a new dataframe to help summarize the findings
ind <- c(1: dim(malignant)[1])
d_sum <- data.frame(ind, de.tgw$table, "twd"=d$tagwise.dispersion)

#plot the tagwise dispersions
hist(d_sum$twd, breaks=20, main = "Tagwise dispersions", xlab = "twd")

#I want the fdr adjusting the pvalue \
#using Benjamini & Hochberg (1995) ("BH" or its alias "fdr") method
d_sum$PValue_fdr <- p.adjust(method="fdr",p=d_sum$PValue)

d_sum[d_sum$Sig_1pct_fdr==TRUE, ]
d_sum[d_sum$Sig_1pct==TRUE, ]

#to write out this useful data frame
write.table(d_sum, file="results_summary.tsv", quote=F)
sessionInfo()
