# exercise

Next Generation Sequencing systems produce count data which can be used for differential testing, with longitudinal studies involving paired samples testing.

Beta-binomial tests have therefore emerged as the methods of choice for the analysis of NGS and also proteomics data involving count data.

In this exercise the user is presented with two groups of samples from the same individual. Sampling was performed in August and in December. For each group sampling was performed 4 times (replicates) (August_1,August_2,August_3,August_4 ; December_1,December_2,December_3,December_4). 

Total counts(t) for a gene as well as counts for a malignant isoform (m) of the same gene are presented in the file raw_data.tsv - (August_1m,August_2m,August_3m,August_4m ; August_1t,August_2t,August_3t,August_4t ; December_1m,December_2m,December_3m,December_4m ; December_1t,December_2t,December_3t,December_4t). Missing values in the raw data correspond to samples that got lost by the logistics partner.

Each row corresponds to a different individual.

The total counts (t) attributed to a gene can for the purpose of this exercise be assumed as the sum of the counts of each of its isoforms i.e., gene (t) = isoform_1 + isoform_2 + isoform_3 + isoform_i. In this exercise: gene(t) = isoform_malignant(m) + isoform_normal.

The objective of this exercise is to identify individuals for which the levels of the malignant isoform (m) significantly changed from August to December independently of the changes in the gene(t).

The solution should be implemented in Python and/or R.
