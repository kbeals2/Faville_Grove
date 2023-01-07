#### ITS AMPLICON PROCESSING FOR FAVILLE GROVE PROJECT ####
# POST PRIMER REMOVAL STEPS (to be completed using supercomputer) #
# Soil samples collected from Faville Prairie (remnant) and Tillotson Prairie (restoration) in summer 2022
# Used AMF specific primers:  WANDA Forward and AML2 Reverse primers

(.packages())

# Installing DADA2 via devtools (necessary when using with AWS RStudio AMI)
install.packages("devtools")
library(devtools)

devtools::install_github("benjjneb/dada2", ref = "v1.18") # note that the AMI created by Louis Aslett is R version 4.0, and the more recent versions of dada2 only work with R 4.2. So we install an earlier version of dada2 here when connecting through AWS AMI

library(dada2)

# Link Dropbox where sequences are stored
library("RStudioAMI")
linkDropbox()
excludeSyncDropbox("*")
includeSyncDropbox("SDSU")


#### 6) Set working directory ####
path <- "/home/rstudio/Dropbox/SDSU/Faville_Grove_microbial_sequences_2022/AMF_seqs_2022/cutadapt2"
list.files(path)

setwd("/home/rstudio/Dropbox/SDSU")


#### 7) Combine all forward and reverse reads. Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_cut <- sort(list.files(path, pattern ="_R1.fastq.gz", full.names = TRUE))
REV_cut <- sort(list.files(path, pattern ="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
(sample.names <- sapply(strsplit(basename(FWD_cut), "_"), `[`, 2)) #sample ID occurs before the first _ in the fastq file


#### 8) Inspect read quality profiles ####
plotQualityProfile(FWD_cut[54:57])
# These sequences are 300 bp in length
# Quality drops below PHRED score of 30 around 200 bases

plotQualityProfile(REV_cut[54:57])
# Quality drops below PHRED score of 30 around 150 bases


#### 9) Filter and Trim ####
# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Since the primers are already removed, we don't truncate here. Instead, we enforce a minimum sequence length (minLen) to get rid of spurious very low-length sequences. 

FWD_filt <- file.path(path, "filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
REV_filt <- file.path(path, "filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(FWD_filt) <- sample.names
names(REV_filt) <- sample.names

out <- filterAndTrim(FWD_cut, FWD_filt, REV_cut, REV_filt, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Ran into a problem when I checked the first 10 samples; only 12% of reads retained after filtering. 

# Try filtering with just the FWD reads b/c the REV reads are lower quality
out2 <- filterAndTrim(FWD_cut, FWD_filt, maxN = 0, maxEE = c(2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Now 59% of reads are retained


#### 10) Learn the error rates ####
# Use the filtered reads to learn the sequence error rates and correct for these in the later steps of the pipeline
err_FWD <- learnErrors(FWD_filt, multithread = TRUE) # took less than 2 min

# Visualize the estimated error rates as a sanity check
# red line = expected error rate based on quality score (note that this decreases as the quality increases)
# black line = estimated error rate after convergence of the machine-learning algorithm
# black dots = observed error frequency in our samples
# black dots should align with black line
plotErrors(err_FWD, nominalQ = TRUE)


#### 11) Sample inference ####
# This algorithm will tell us how many unique sequences are in each sample after controlling for sequencing errors. Can also set pool = TRUE to pool across all samples to detect low abundance variants.
dada_FWD <- dada(FWD_filt, err = err_FWD, multithread = TRUE, pool = "pseudo")


#### 12) Construct a sequence table ####
AMF_seq_table <- makeSequenceTable(dada_FWD)
dim(AMF_seq_table) 
# (number of samples, number of amplicon sequence variants ASVs) (found 6,277 unique sequences across the 57 samples)


#### 13) Remove chimeras ###
AMF_seq_table_nochim <- removeBimeraDenovo(AMF_seq_table, method = "consensus", multithread = TRUE, verbose = TRUE)
# found 1,020 chimeras out of 6,277 unique sequences (16% identified as chimeric)

sum(AMF_seq_table_nochim)/sum(AMF_seq_table)
# retained 98.9% of sequence abundance


#### 14) Track reads through the pipeline to verify everything worked as expected ####
getN <- function(x) sum(getUniques(x))
(summary_table <- data.frame(row.names = sample.names, 
                             input = out[, 1],
                             filtered = out[, 2], 
                             denoised_FWD = sapply(dada_FWD, getN),
                             non_chim = rowSums(AMF_seq_table_nochim), 
                             final_perc_reads_retained = round(rowSums(AMF_seq_table_nochim)/out[, 1]*100, 1)))
# Looks ok. On average, retained __% of reads. 


#### 15) Make summary tables that can be used to assign taxonomy

# giving seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(AMF_seq_table_nochim)
asv_headers <- vector(dim(AMF_seq_table_nochim)[2], mode = "character")

for (i in 1:dim(AMF_seq_table_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# make a fasta file of the final ASV seqs (once this is made, upload to rdp.cme.msu.edu/classifier/classifier.jsp)
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/rstudio/Dropbox/SDSU/FG_2022_AMF_ASV.fa")

# make count table
asv_table <- t(AMF_seq_table_nochim)
row.names(asv_table) <- sub(">", "", asv_headers)
write.table(asv_table, "FG_2022_AMF_ASV_counts.tsv", sep = "\t", quote = F, col.names = NA)  