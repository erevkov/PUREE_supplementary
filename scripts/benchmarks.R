library(plyr)
library(tidyverse)

setwd("/home/ubuntu/purity-ml/")
### Read purity estimates of other methods
## DeMixT and LinSeed use counts, CIBERSORTx prefers TPM

# table with protein coding genes - their ENSEMBL IDs and HGNC symbols
protein.coding.genes.names <- read.csv("./data/proteinc_gene_name_conversion_table.csv", row.names = 1)
test.samples.names <- as.vector(read.csv("./data/test_samples_names.tsv")[, 1])
train.samples.names <- as.vector(read.csv("./data/train_samples_names.tsv")[, 1])

samples.names <- train.samples.names

### data loading
### TCGA TPM data
log2.tpm_pancan <- read.csv("./data/tcga_tpm_pancan_60k_10535.tsv", sep='\t', 
                            row.names = 1, check.names = FALSE)
rownames(log2.tpm_pancan) <- sapply(strsplit(rownames(log2.tpm_pancan),"\\."), `[`, 1) # cutting version postfixes of ensemble ids

## only for TPM values
tumor.log2.tpm <- log2.tpm_pancan[, samples.names]
tumor.log2.tpm <- tumor.log2.tpm[rownames(tumor.log2.tpm) %in% protein.coding.genes.names$ensembl_id, ]
tumor.tpm <- 2**tumor.log2.tpm - 0.001
hgnc_tpm_gene_names <- mapvalues(rownames(tumor.tpm),
                                 from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.tpm), ]$ensembl_id),
                                 to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.tpm), ]$hgnc_sym))


### TCGA counts data (TOIL RSEM expected_count)
log2.counts <- read.csv("./data/tcga_gene_expected_count.tsv", sep = '\t', 
                        row.names = 1, check.names = FALSE)

#### TCGA: choose the code below depending on which data you have loaded above
# data frame to matrix conversion + formatting
log2.counts = data.matrix(log2.counts) 
rownames(log2.counts) <- sapply(strsplit(rownames(log2.counts),"\\."), `[`, 1) # cutting version postfixes of ensemble ids 
# leave only protein-coding genes 
log2.counts <- log2.counts[rownames(log2.counts) %in% protein.coding.genes.names$ensembl_id, ]
# move to linear scale, log2(x+1)
counts <- 2**log2.counts - 1
# leave only test samples
normal.counts <- (counts[, grep('.11$',colnames(counts)), drop=FALSE]) # get all possible normal samples (we're only using them for DeMixT)
common.names <- intersect(colnames(counts), samples.names) # counts doesn't have "TCGA-SW-A7EB-01"? 
tumor.counts <- (counts[, common.names, drop=FALSE]) 
# names of genes in HGNC nomenclature
hgnc_counts_gene_names <- mapvalues(rownames(tumor.counts), 
                                    from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.counts), ]$ensembl_id), 
                                    to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.counts), ]$hgnc_sym))

### Chen et al Lung data 
## fpkm
# chen.tumor.fpkm <- read.csv("./data/Chen_et_al/RSEM_FPKM_tumor.tsv", sep = '\t', 
#                            row.names = 1, check.names = FALSE)
# chen.tumor.fpkm <- chen.tumor.fpkm[rownames(chen.tumor.fpkm) %in% protein.coding.genes.names$ensembl_id, ]
chen.tumor.tpm <- read.csv("./data/Chen_et_al/TPM_tumor.tsv", sep = '\t', 
                          row.names = 1, check.names = FALSE)
chen.tumor.tpm <- chen.tumor.tpm[rownames(chen.tumor.tpm) %in% protein.coding.genes.names$ensembl_id, ]
hgnc_gsk_gene_names <- mapvalues(rownames(chen.tumor.tpm),
                                 from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(chen.tumor.tpm), ]$ensembl_id),
                                 to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(chen.tumor.tpm), ]$hgnc_sym))

## counts
# chen.tumor.counts <- read.csv("./data/Chen_et_al/RSEM_counts_tumor.tsv", sep = '\t', 
#                                   row.names = 1, check.names = FALSE)
chen.tumor.counts <- read.csv("./data/Chen_et_al/salmon_counts_tumor.tsv", sep = '\t',
                             row.names = 1, check.names = FALSE)
chen.tumor.counts <- chen.tumor.counts[rownames(chen.tumor.counts) %in% protein.coding.genes.names$ensembl_id, ]
chua.normal.counts <- read.csv("./data/Chen_et_al/RSEM_counts_normal.tsv", sep = '\t', 
                              row.names = 1, check.names = FALSE)
chua.normal.counts <- chua.normal.counts[rownames(chua.normal.counts) %in% protein.coding.genes.names$ensembl_id, ]
# adjust the purity computation lines below accordingly

### Chua et al. - Lung 
chua.tumor.fpkm <- read.csv('./data/Chua_et_al/FPKM_tumor.csv', row.names = 1, check.names = FALSE) # rows already in HGNC ids
chua.tumor.fpkm <- chua.tumor.fpkm[rownames(chua.tumor.fpkm) %in% protein.coding.genes.names$hgnc_sym, ] 
chua.tumor.fpkm <- chua.tumor.fpkm ** 2
hgnc_tki_gene_names <- mapvalues(rownames(chua.tumor.fpkm),
                                 from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(chua.tumor.fpkm), ]$ensembl_id),
                                 to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(chua.tumor.fpkm), ]$hgnc_sym))

### Juanito et al. # linear space by default
juanito.tumor.tpm <- read.csv('./data/Joanito_et_al/TPM.tsv', sep = '\t', row.names = 2, check.names = FALSE)[,-c(1,2)] # rows already in HGNC ids
juanito.tumor.tpm <- juanito.tumor.tpm[rownames(juanito.tumor.tpm) %in% protein.coding.genes.names$hgnc_sym, ] 
hgnc_nant_gene_names <- mapvalues(rownames(juanito.tumor.tpm),
                                  from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(juanito.tumor.tpm), ]$ensembl_id),
                                  to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(juanito.tumor.tpm), ]$hgnc_sym))

### TCGA (non-PUREE samples)
# TGCT (FPKM)
tcga_tgct.tumor.log2.fpkm <- read.csv('./data/TCGA_new/TCGA-TGCT_htseq_fpkm_formatted.tsv', sep='\t',
                                 row.names = 1, check.names = FALSE)
tcga_tgct.tumor.fpkm <- t(2**tcga_tgct.tumor.log2.fpkm - 1)
tcga_tgct.tumor.fpkm <- tcga_tgct.tumor.fpkm[rownames(tcga_tgct.tumor.fpkm) %in% protein.coding.genes.names$ensembl_id, ]
# TGCT (Counts)
tcga_tgct.tumor.log2.counts <- read.csv('./data/TCGA_new/TCGA-TGCT_htseq_counts_formatted.tsv', sep='\t',
                                      row.names = 1, check.names = FALSE)
tcga_tgct.tumor.counts <- t(2**tcga_tgct.tumor.log2.counts - 1)
tcga_tgct.tumor.counts <- tcga_tgct.tumor.counts[rownames(tcga_tgct.tumor.counts) %in% protein.coding.genes.names$ensembl_id, ]
# CRC (FPKM)
tcga_crc.tumor.log2.fpkm <- read.csv('./data/TCGA_new/TCGA-CRC_new_htseq_fpkm_formatted.tsv', sep='\t',
                                      row.names = 1, check.names = FALSE) 
tcga_crc.tumor.fpkm <- t(2**tcga_crc.tumor.log2.fpkm - 1)
tcga_crc.tumor.fpkm <- tcga_crc.tumor.fpkm[rownames(tcga_crc.tumor.fpkm) %in% protein.coding.genes.names$ensembl_id, ]
# CRC (Counts)
tcga_crc.tumor.log2.counts <- read.csv('./data/TCGA_new/TCGA-CRC_new_htseq_counts_formatted.tsv', sep='\t',
                                        row.names = 1, check.names = FALSE)
tcga_crc.tumor.counts <- t(2**tcga_crc.tumor.log2.counts - 1)
tcga_crc.tumor.counts <- tcga_crc.tumor.counts[rownames(tcga_crc.tumor.counts) %in% protein.coding.genes.names$ensembl_id, ]
## normals
tcga_crc.normal.log2.counts <- read.csv('./data/TCGA_new/TCGA-CRC_new_normals_htseq_counts_formatted.tsv', sep='\t',
                                        row.names = 1, check.names = FALSE)
tcga_crc.normal.counts <- t(2**tcga_crc.normal.log2.counts - 1)
tcga_crc.normal.counts <- tcga_crc.normal.counts[rownames(tcga_crc.normal.counts) %in% protein.coding.genes.names$ensembl_id, ]
# UCEC (FPKM)
tcga_ucec.tumor.log2.fpkm <- read.csv('./data/TCGA_new/TCGA-UCEC_new_htseq_fpkm_formatted.tsv', sep='\t',
                                     row.names = 1, check.names = FALSE) 
tcga_ucec.tumor.fpkm <- t(2**tcga_ucec.tumor.log2.fpkm - 1)
tcga_ucec.tumor.fpkm <- tcga_ucec.tumor.fpkm[rownames(tcga_ucec.tumor.fpkm) %in% protein.coding.genes.names$ensembl_id, ]
# UCEC (Counts)
tcga_ucec.tumor.log2.counts <- read.csv('./data/TCGA_new/TCGA-UCEC_new_htseq_counts_formatted.tsv', sep='\t',
                                      row.names = 1, check.names = FALSE) 
tcga_ucec.tumor.counts <- t(2**tcga_ucec.tumor.log2.counts - 1)
tcga_ucec.tumor.counts <- tcga_ucec.tumor.counts[rownames(tcga_ucec.tumor.counts) %in% protein.coding.genes.names$ensembl_id, ]
## normals
tcga_ucec.normal.log2.counts <- read.csv('./data/TCGA_new/TCGA-UCEC_new_normals_htseq_counts_formatted.tsv', sep='\t',
                                        row.names = 1, check.names = FALSE) 
tcga_ucec.normal.counts <- t(2**tcga_ucec.normal.log2.counts - 1)
tcga_ucec.normal.counts <- tcga_ucec.normal.counts[rownames(tcga_ucec.normal.counts) %in% protein.coding.genes.names$ensembl_id, ]
# TGCT (FPKM)
tcga_pcpg.tumor.log2.fpkm <- read.csv('./data/TCGA_new/TCGA-PCPG_htseq_fpkm_formatted.tsv', sep='\t',
                                      row.names = 1, check.names = FALSE)
tcga_pcpg.tumor.fpkm <- t(2**tcga_pcpg.tumor.log2.fpkm - 1)
tcga_pcpg.tumor.fpkm <- tcga_pcpg.tumor.fpkm[rownames(tcga_pcpg.tumor.fpkm) %in% protein.coding.genes.names$ensembl_id, ]
# TGCT (Counts)
tcga_pcpg.tumor.log2.counts <- read.csv('./data/TCGA_new/TCGA-PCPG_htseq_counts_formatted.tsv', sep='\t',
                                        row.names = 1, check.names = FALSE)
tcga_pcpg.tumor.counts <- t(2**tcga_pcpg.tumor.log2.counts - 1)
tcga_pcpg.tumor.counts <- tcga_pcpg.tumor.counts[rownames(tcga_pcpg.tumor.counts) %in% protein.coding.genes.names$ensembl_id, ]

###########
### CibersortX 
# uses its own online interface, requires non-log data and gene names in HGNC symbol nomenclature
# signature HNSCC matrix is found in CibersortX's paper supplementary 2e (named as 2d)
# signature NSCLC matrix is found in CibersortX's paper supplementary 2l 
data.cibersort <- tcga_pcpg.tumor.fpkm
hgnc_gene_names <- mapvalues(rownames(data.cibersort),
                                  from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(data.cibersort), ]$ensembl_id),
                                  to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(data.cibersort), ]$hgnc_sym))

data.cibersort.df <- rbind(setNames(data.frame(t(c("gene", colnames(data.cibersort))), stringsAsFactors = FALSE), c("gene", colnames(data.cibersort))), 
                           setNames(data.frame(cbind(hgnc_gene_names, data.cibersort), check.names = FALSE, stringsAsFactors = FALSE), c("gene", colnames(data.cibersort))))
# saving results
# write.table(data.cibersort.df, file = "./data/set_data/CibersortX_test_set_samples_tpm.tsv", sep='\t',
# quote=FALSE, col.names=FALSE, row.names=FALSE)

##########
# Cibersort's purities
# cibersort.results <- read.csv("./results/CIBERSORTx_test_results.csv", row.names=1)

cibersort.purity <- select(cibersort.results, EPCAM)
# write.csv(cibersort.purity, "./results/CIBERSORTx_test_purities.csv")

###########
### LinSeed
# devtools::install_github("ctlab/linseed")
library(linseed)

# normalization recommended by the package authors: linear space + normalized column-wise so they the columns (samples) have the same sum
data.linseed <- tcga_pcpg.tumor.counts
hgnc_gene_names <- mapvalues(rownames(data.linseed),
                             from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(data.linseed), ]$ensembl_id),
                             to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(data.linseed), ]$hgnc_sym))
expData <- data.linseed
rownames(expData) <- hgnc_gene_names
expData <- expData[!grepl("RPL|RPS", rownames(expData)), ] # removing RPL|RPS genes
expData <- data.matrix(expData)
expData[is.nan(expData)] <- 0
# expData[is.na(expData)] <- 0
totalSum <- mean(colSums(expData))
expData <- apply(expData, 2, function(x) x * totalSum / sum(x))

lo <- linseed::LinseedObject$new(expData, topGenes=10000)

run_linseed <- function(lo) {
  # Colinearity networks
  # the dataset has to be in lin-scale and normalized
  lo$calculatePairwiseLinearity()
  lo$calculateSpearmanCorrelation()
  lo$calculateSignificanceLevel(100)
  # lo$significancePlot(0.01)
  
  lo$filterDatasetByPval(0.01)
  
  lo$setCellTypeNumber(2)
  lo$project("filtered")
  # lo$projectionPlot(color="filtered")
  lo$smartSearchCorners(dataset="filtered", error="norm")
  # deconvolution
  lo$deconvolveByEndpoints()
  # lo$selectGenes(100)
  # lo$tsnePlot()
  return(lo)
}
lo <- run_linseed(lo)
linseed.purity <- lo$proportions[2, ] # 2nd comp is purity judging by the values / good correlation

write.csv(t(lo$proportions), "./results/LinSeed_TCGA_PCPG_proportions.csv")

#### TIME/MEMORY BENCHMARK ####
linseed_benchmark <- bench::mark({run_linseed(lo)}, iterations=1)
run_time <- as.numeric(linseed_benchmark$total_time) # sec
total_mem_alloc <- as.numeric(linseed_benchmark$mem_alloc)/(1024^2) # mb
max_mem <- max(linseed_benchmark$memory[[1]]$bytes, na.rm=TRUE)/(1024^2) # mb
data.frame(run_time, total_mem_alloc, max_mem, cohort='TCGA_PCPG', method='LinSeed')  %>%
  write.csv(., './results/time_memory_benchmark/LinSeed_prof_benchmark_TCGA_PCPG.csv')
# TCGA test set: 28455.8 MB, 192980 ms, IO:
#### TIME/MEMORY BENCHMARK END ####

# saving results
# write.csv(linseed.purity, "./results/LinSeed_purities.csv")
###########

###########
### DeMixT-normalization
# quartile normalization provided by the authors of DeMixT in their email
# BiocManager::install("DSS")
library(DSS)

quart_normalize_for_demixt <- function(normal.counts, tumor.counts) {
  
  Count.matrix <- cbind(normal.counts, tumor.counts)
  newt <- Count.matrix
  colnames(newt)=NULL
  rownames(newt)=NULL
  
  normal.id <- colnames(normal.counts)
  tumor.id <- colnames(tumor.counts)
  
  designs=c(rep("1", length(normal.id)), rep("0", length(tumor.id)))
  seqData=newSeqCountSet(as.matrix(newt), designs)
  
  # Quartile normalization/total/median (median is provided by the authors, seems to work fine)
  seqData=estNormFactors(seqData, "quantile")
  k3=seqData@normalizationFactor
  mk3=median(k3)
  k3=k3/mk3
  
  temp <- newt
  
  for(i in 1:ncol(newt)){
    temp[,i] = temp[,i]/k3[i]
  }
  Count.matrix.normalized <- temp
  colnames(Count.matrix.normalized) <- colnames(Count.matrix)
  rownames(Count.matrix.normalized) <- rownames(Count.matrix)
  
  normal.counts.quartnorm <- Count.matrix.normalized[, colnames(normal.counts), drop=FALSE]
  tumor.counts.quartnorm <- Count.matrix.normalized[, colnames(tumor.counts), drop=FALSE]
  
  return(list("normal" = normal.counts.quartnorm, 
              "tumor" = tumor.counts.quartnorm))
}

quartnorm_counts_list <- quart_normalize_for_demixt(tcga_ucec.normal.counts, tcga_ucec.tumor.counts)
normal.counts.quartnorm <- quartnorm_counts_list$normal
tumor.counts.quartnorm <- quartnorm_counts_list$tumor

# only leave the genes where the values are >= 1
tumor.genes <- rownames(tumor.counts.quartnorm[rowSums(tumor.counts.quartnorm < 1) == 0, ])
normal.genes <- rownames(normal.counts.quartnorm[rowSums(normal.counts.quartnorm < 1) == 0, ])
common.demixt.genes <- intersect(tumor.genes, normal.genes)

### DeMixT-run
# Note: author's used quartile-normalized counts data, although they say 
#   that theoretically DeMixT should accept any input expression as long as the distribution follows
#   log2-normal distribution
#devtools::install_github("wwylab/DeMixTallmaterials/DeMixT_0.2")
# OR
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# OR (recommended on github)
# BiocManager::install("DeMixT")
library(DeMixT)

res <- DeMixT(data.Y = SummarizedExperiment(tumor.counts.quartnorm[common.demixt.genes, ]), 
              data.N1 = SummarizedExperiment(normal.counts.quartnorm[common.demixt.genes, ])
)
demixt.stromal <- res$pi[1,]
demixt.purity <- 1 - demixt.stromal
# saving results
# write.csv(demixt.purity, "./results/DeMixT_test_purities.csv")
### TIME MEMORY BENCHMARK
demixt_benchmark <- bench::mark({DeMixT(data.Y = SummarizedExperiment(tumor.counts.quartnorm[common.demixt.genes, ]), 
                                        data.N1 = SummarizedExperiment(normal.counts.quartnorm[common.demixt.genes, ])
)}, iterations=1)
run_time <- as.numeric(demixt_benchmark$total_time) # sec
total_mem_alloc <- as.numeric(demixt_benchmark$mem_alloc)/(1024^2) # mb
max_mem <- max(demixt_benchmark$memory[[1]]$bytes, na.rm=TRUE)/(1024^2) # mb
data.frame(run_time, total_mem_alloc, max_mem, cohort='TCGA_UCEC_new', method='DeMixT')  %>%
  write.csv(., './results/time_memory_benchmark/DeMixT_prof_benchmark_TCGA_UCEC_new.csv')
### TIME MEMORY BENCHMARK END
###########


########### EPIC
# devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
library(EPIC)
gene_expr.epic <- tcga_pcpg.tumor.fpkm
hgnc_gene_names_epic <- mapvalues(rownames(gene_expr.epic),
                                  from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(gene_expr.epic), ]$ensembl_id),
                                  to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(gene_expr.epic), ]$hgnc_sym))
rownames(gene_expr.epic) <- hgnc_gene_names_epic
# purity calculation
out <- EPIC(bulk = gene_expr.epic)
epic.purity <- out$cellFractions[,'otherCells'] # otherCells~cancer cells in tumor samples, though might also include stromal, epithelial or others
# write.csv(epic.purity, "./results/EPIC_test_purities.csv")

#### TIME/MEMORY BENCHMARK ####
epic_benchmark <- bench::mark(EPIC(bulk = gene_expr.epic), iterations=1)
run_time <- as.numeric(epic_benchmark$total_time) # sec
total_mem_alloc <- as.numeric(epic_benchmark$mem_alloc)/(1024^2) # mb
max_mem <- max(epic_benchmark$memory[[1]]$bytes, na.rm=TRUE)/(1024^2) # mb
data.frame(run_time, total_mem_alloc, max_mem, cohort='TCGA_PCPG', method='EPIC')  %>%
  write.csv(., './results/time_memory_benchmark/EPIC_prof_benchmark_TCGA_PCPG.csv')
#### TIME/MEMORY BENCHMARK END ####

##### ESTIMATE
# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)

# setwd("/home/ubuntu/purity-ml/ESTIMATE")

estimate.rnaseq.data <- tcga_pcpg.tumor.fpkm
hgnc_gene_names <- mapvalues(rownames(estimate.rnaseq.data),
                                  from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(estimate.rnaseq.data), ]$ensembl_id),
                                  to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(estimate.rnaseq.data), ]$hgnc_sym))
rownames(estimate.rnaseq.data) <- hgnc_gene_names

input.rnaseq.estimate <- "./ESTIMATE/RNA-seq_gene_expression_hgnc.tsv"
gct.rnaseq <- "./ESTIMATE/RNA-seq_gene_expression.gct"
gct.rnaseq.filtered <- "./ESTIMATE/RNA-seq_gene_expression_filtered.gct"
scores.path <- "./ESTIMATE/RNA-seq_gene_expression_ESTIMATE_scores.tsv"

# preparing data for ESTIMATE
write.table(estimate.rnaseq.data, file=input.rnaseq.estimate, row.names = TRUE, sep = "\t")

run_estimate <- function() {
  outputGCT_fixed(input.rnaseq.estimate, gct.rnaseq) # BUGGY in original version: quote="" when reading the table makes the samples names to bug out later
  filterCommonGenes_fixed(gct.rnaseq, gct.rnaseq.filtered, id="GeneSymbol") # BUGGY in original version because of outputGCT 
  estimateScore(gct.rnaseq.filtered, scores.path, platform="affymetrix")
}

run_estimate()

estimate.scores <- read.csv2("./ESTIMATE/RNA-seq_gene_expression_ESTIMATE_scores.tsv", sep='\t')
# creating df with purities
sample_ids <- gsub("\\.", "-", unlist(as.vector(estimate.scores[2, 3:length(estimate.scores)])))
estimate.purity <- estimate.scores[6, 3:length(estimate.scores)]
colnames(estimate.purity) <- sample_ids 
estimate.purity <- data.frame(t(estimate.purity))
colnames(estimate.purity) <- 'ESTIMATE'

write.csv(estimate.purity, "./results/ESTIMATE_TCGA_PCPG_purities.csv")

# after the fix the scores look ok...
### TIME MEMORY BENCHMARK
estimate_benchmark <- bench::mark(run_estimate())
run_time <- as.numeric(estimate_benchmark$total_time) # sec
total_mem_alloc <- as.numeric(estimate_benchmark$mem_alloc)/(1024^2) # mb
max_mem <- max(estimate_benchmark$memory[[1]]$bytes, na.rm=TRUE)/(1024^2) # mb
data.frame(run_time, total_mem_alloc, max_mem, cohort='TCGA_PCPG', method='ESTIMATE')  %>%
  write.csv(., './results/time_memory_benchmark/ESTIMATE_prof_benchmark_TCGA_PCPG.csv')
### TIME MEMORY BENCHMARK END

### FIXED ESTIMATE FUNCTIONS
# FIXED outputGCT:
# input.f <- input.rnaseq
# output.f <- gct.rnaseq

outputGCT_fixed <- function(input.f, output.f){
  if (is.data.frame(input.f) == TRUE) {
    exp.data <- input.f
  }
  else {
    # exp.data <- read.table(input.f, header = TRUE, row.names = 1, 
    #                        sep = "\t", quote = "")
    exp.data <- read.table(input.f, header = TRUE, row.names = 1, sep = "\t")
  }
  exp.data1 <- data.frame(NAME = rownames(exp.data), Description = rownames(exp.data), 
                          exp.data, stringsAsFactors=FALSE)
  column1 <- colnames(exp.data1)
  column1[1] <- "NAME"
  column1[2] <- "Description"
  exp.data1$NAME <- factor(exp.data1$NAME)
  exp.data1$Description <- factor(exp.data1$Description)
  levels(exp.data1[, 1]) <- c(levels(exp.data1[, 1]), "NAME")
  levels(exp.data1[, 2]) <- c(levels(exp.data1[, 2]), "Description")
  exp.data2 <- rbind(column1, exp.data1)
  exp.data2 <- rbind(column1, data.frame(exp.data1, stringsAsFactors = FALSE))
  row1 <- rep("", length(1:ncol(exp.data)))
  row1_2 <- data.frame(row1, row1)
  row1_2 <- t(row1_2)
  No_gene <- nrow(exp.data1)
  No_sample <- (ncol(exp.data1) - 2)
  GCT <- matrix(c("#1.2", No_gene, "", No_sample), nrow = 2, 
                ncol = 2)
  gct <- cbind(GCT, row1_2)
  colnames(gct) <- colnames(exp.data2)
  tmp <- rbind(gct, exp.data2)
  write.table(tmp, output.f, sep = "\t", row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  invisible(NULL)
}

input.f <- gct.rnaseq
output.f <- gct.rnaseq.filtered
id <- "GeneSymbol"

filterCommonGenes_fixed <- function(input.f, output.f, id = c("GeneSymbol", "EntrezID")){
  stopifnot((is.character(input.f) && length(input.f) == 1 && 
               nzchar(input.f)) || (inherits(input.f, "connection") && 
                                      isOpen(input.f, "r")))
  stopifnot((is.character(output.f) && length(output.f) == 
               1 && nzchar(output.f)))
  id <- match.arg(id)
  input.df <- read.table(input.f, header = TRUE, row.names = 1, 
                         sep = "\t", quote = "", stringsAsFactors = FALSE)
  merged.df <- merge(common_genes, input.df, by.x = id, by.y = "row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  merged.df <- merged.df[, -1:-ncol(common_genes)]
  colnames(merged.df) <- input.df[1, ] # FIX: or it won't write the samples ids
  merged.df <- merged.df[, -1] # FIX: remove unnecessary description column
  print(sprintf("Merged dataset includes %d genes (%d mismatched).", 
                nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
  outputGCT_fixed(merged.df, output.f)
}

function (input.ds, output.ds, platform = c("affymetrix", "agilent", 
                                            "illumina")) 
{
  stopifnot(is.character(input.ds) && length(input.ds) == 1 && 
              nzchar(input.ds))
  stopifnot(is.character(output.ds) && length(output.ds) == 
              1 && nzchar(output.ds))
  platform <- match.arg(platform)
  ds <- read.delim(input.ds, header = TRUE, sep = "\t", skip = 2, 
                   row.names = 1, blank.lines.skip = TRUE, as.is = TRUE, 
                   na.strings = "")
  descs <- ds[, 1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  dataset <- list(ds = ds, row.names = row.names, descs = descs, 
                  names = names)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/")
  s <- length(temp[[1]])
  input.file.name <- temp[[1]][s]
  temp <- strsplit(input.file.name, split = ".gct")
  input.file.prefix <- temp[[1]][1]
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  gs <- as.matrix(SI_geneset[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(SI_geneset)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", 
                length(gene.overlap)))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    }
    else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG/Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector/sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        }
        else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES = ES, arg.ES = arg.ES, 
                                RES = RES, indicator = TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  estimate.score <- apply(score.data, 2, sum)
  if (platform != "affymetrix") {
    score.data <- rbind(score.data, estimate.score)
    rownames(score.data) <- c("StromalScore", "ImmuneScore", 
                              "ESTIMATEScore")
  }
  else {
    convert_row_estimate_score_to_tumor_purity <- function(x) {
      stopifnot(is.numeric(x))
      cos(0.6049872018 + 0.0001467884 * x)
    }
    est.new <- NULL
    for (i in 1:length(estimate.score)) {
      est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
      est.new <- rbind(est.new, est_i)
      if (est_i >= 0) {
        next
      }
      else {
        message(paste(sample.names[i], ": out of bounds", 
                      sep = ""))
      }
    }
    colnames(est.new) <- c("TumorPurity")
    estimate.t1 <- cbind(estimate.score, est.new)
    x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 
      0
    estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
    score.data <- rbind(score.data, t(estimate.t1))
    rownames(score.data) <- c("StromalScore", "ImmuneScore", 
                              "ESTIMATEScore", "TumorPurity")
  }
  outputGCT(score.data, output.ds)
}

##### DeconRNASeq
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DeconRNASeq")
library(DeconRNASeq)

gene_expr.deconrnaseq <- as.data.frame(tcga_pcpg.tumor.fpkm)
hgnc_gene_names_deconrnaseq <- mapvalues(rownames(gene_expr.deconrnaseq),
                                         from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(gene_expr.deconrnaseq), ]$ensembl_id),
                                         to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(gene_expr.deconrnaseq), ]$hgnc_sym))
rownames(gene_expr.deconrnaseq) <- hgnc_gene_names_deconrnaseq

signature <- read.csv2('./data/CIBERSORTx_NSCLC_signature.tsv', sep='\t', row.names=1)
sig <- mutate_all(signature, function(x) as.numeric(as.character(x)))
rownames(sig) <- rownames(signature)

res <- DeconRNASeq(datasets=gene_expr.deconrnaseq, signatures=sig)
deconrnaseq.purity <- data.frame(res[[1]][,'EPCAM'], row.names=colnames(gene_expr.deconrnaseq))
colnames(deconrnaseq.purity) <- 'DeconRNASeq'
write.csv(deconrnaseq.purity, "./results/DeconRNASeq_TCGA_PCPG_purities.csv")

#### TIME/MEMORY BENCHMARK ####
deconrnaseq_benchmark <- bench::mark(DeconRNASeq(datasets=gene_expr.deconrnaseq, signatures=sig))
run_time <- as.numeric(deconrnaseq_benchmark$total_time) # sec
total_mem_alloc <- as.numeric(deconrnaseq_benchmark$mem_alloc)/(1024^2) # mb
max_mem <- max(deconrnaseq_benchmark$memory[[1]]$bytes, na.rm=TRUE)/(1024^2) # mb
data.frame(run_time, total_mem_alloc, max_mem, cohort='TCGA_PCPG', method='DeconRNASeq')  %>%
  write.csv(., './results/time_memory_benchmark/DeconRNASeq_prof_benchmark_TCGA_PCPG.csv')
#### TIME/MEMORY BENCHMARK END ####