#### 16S AMPLICON PROCESSING FOR FAVILLE GROVE PROJECT ####
# Soil samples collected from Faville Prairie (remnant) and Tillotson Prairie (restoration) in summer 2022
# Used the Earth Microbiome Project primers: 515F/806R


# Installing via Bioconductor
(.packages())

install.packages("BiocManager")
library(BiocManager)

BiocManager::install("Biostrings")
library(Biostrings)

BiocManager::install("ShortRead")
library(ShortRead)

library(dada2)
packageVersion("dada2")


# Installing via devtools (necessary when using with AWS RStudio AMI)
install.packages("devtools")
library(devtools)

devtools::install_github("benjjneb/dada2", ref = "v1.18") # note that the AMI created by Louis Aslett is R version 4.0, and the more recent versions of dada2 only work with R 4.2. So we install an earlier version of dada2 here when connecting through AWS AMI

library(dada2)

# Link Dropbox where sequences are stored
library("RStudioAMI")
linkDropbox()
excludeSyncDropbox("*")
includeSyncDropbox("SDSU")


#### 1) Set working directory ####
setwd("/home/rstudio/Dropbox/SDSU")

path <- "/home/rstudio/Dropbox/SDSU/Faville_Grove_microbial_sequences_2022/16S_seqs_2022"
list.files(path)

#### 2) Combine all forward and reverse reads. Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_reads <- sort(list.files(path, pattern ="_R1.fastq.gz", full.names = TRUE))
REV_reads <- sort(list.files(path, pattern ="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
(sample.names <- sapply(strsplit(basename(FWD_reads), "_"), `[`, 2)) #sample ID occurs before the second _ in the fastq file


#### 3) Inspect read quality profiles ####
plotQualityProfile(FWD_reads[1:4])
# Quality doesn't drop below q30; no need to truncate

plotQualityProfile(REV_reads[1:4])
# Quality doesn't drop below q30; no need to truncate


#### 4) Check for primers on the sequences ####
FWD_primer <- "GTGYCAGCMGCCGCGGTAA"
REV_primer <- "GGACTACNVGGGTWTCTAA"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

(FWD_orients <- allOrients(FWD_primer))
(REV_orients <- allOrients(REV_primer))

# The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. We are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
FWD_filtN <- file.path(path, "filtN", basename(FWD_reads)) # Put N-filterd files in filtered/ subdirectory
REV_filtN <- file.path(path, "filtN", basename(REV_reads))
filterAndTrim(FWD_reads, FWD_filtN, REV_reads, REV_filtN, maxN = 0, multithread = TRUE)

# Now count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. It's not necessary to check the presence of the primers on all the samples, so we'll check a few samples randomly.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD_ForwardReads = sapply(FWD_orients, primerHits, fn = FWD_filtN[[1]]), 
      FWD_ReverseReads = sapply(FWD_orients, primerHits, fn = REV_filtN[[1]]), 
      REV_ForwardReads = sapply(REV_orients, primerHits, fn = FWD_filtN[[1]]), 
      REV_ReverseReads = sapply(REV_orients, primerHits, fn = REV_filtN[[1]]))

# As expected, the FWD primer is found in the FWD reads in its forward orientation and on the REV reads in its reverse complement orientation. The REV primer is found in the REV reads in its forward orientation and on the FWD reads in its reverse complement orientation. 

## Before the filter and trim step, need to do some basic math to determine the truncation parameters
# F primer: GTGYCAGCMGCCGCGGTAA (spans base positions 515-533)
# R primer: GGACTACNVGGGTWTCTAA (spans base positions 806-824)
# dada2's default mergePairs command requires a minimum of 12 bases of overlap. To calculate the overlap that we have with our sequences, subtract the target amplicon length from the sum of the truncated FWD & REV reads. 
# The FWD primer starts at base position 515 and the REV primer starts at base position 806. Since the REV primer is read from right to left, subtract the start of the FWD primer from the start of REV primer. The length of the target amplicon is 806 - 515 + 1 = 292 bases.
# Since the FWD and REV reads are great quality, they don't need to be truncated at all, so both reads will stay at 150 bases in length. The FWD and REV reads together sum to 300 bases in length. 300 bases - 292 bases (length of target amplicon) is only 8 bases, which is less than the default parameter of 12 bases of overlap needed to merge FWD and REV reads. We can change the overlap parameter in the mergePairs step later in the pipeline.
# Helpful dada2 support issue: https://github.com/benjjneb/dada2/issues/425

### PREVIOUS RATIONALE (incorrect)
# Target amplicon sequence is 274 bps in length (806-533 + 1)
# We have 300 bps of length with the F & R reads combined (150 bps each), so 300 bps - 274 (target amplicon length) = 26 bps of overlap
# dada2's mergePairs command needs at least 12 bases to be overlapping, but let's say 20 to be conservative. So, 26-20 = 6; 6 divided by 2 = 3 bps could be trimmed from each F & R read
# Ben Callahan recommends that the sum of the truncated FWD & REV reads should be 20 bases longer than the length of the sequenced amplicon. If the reads are truncated to 149 bases each, that sums to 298 bases. 274 bases (length of the amplicon) + 20 = 294.


#### 5) Filter and Trim ####
# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note that the primers are still on the FWD & REV reads, so we remove them with the trimLeft parameter(length of FWD primer = 19 bases, length of REV primer = 19 bases). The reads are in great quality so we don't need to truncate them less than their length of 150 bases.
filt_FWD <- file.path(path, "filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
filt_REV <- file.path(path, "filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(filt_FWD) <- sample.names
names(filt_REV) <- sample.names
out <- filterAndTrim(FWD_reads, filt_FWD, REV_reads, filt_REV, trimLeft = c(19, 19), truncLen = c(150, 150),
                     maxN = 0, maxEE = c(2,5), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE, verbose = TRUE)


#### 6) Inspect read quality profiles after trimming ####
plotQualityProfile(filt_FWD[1:4])
plotQualityProfile(filt_REV[1:4])


#### 7) Learn the error rates ####
# Use the filtered reads to learn the sequence error rates and correct for these in the later steps of the pipeline
err_FWD <- learnErrors(filt_FWD, multithread = TRUE) # took ~ 10 min using AWS AMI
err_REV <- learnErrors(filt_REV, multithread = TRUE) # took ~ 15 min using AWS AMI

# Visualize the estimated error rates as a sanity check
# red line = expected error rate based on quality score (note that this decreases as the quality increases)
# black line = estimated error rate after convergence of the machine-learning algorithm
# black dots = observed error frequency in our samples
# black dots should align with black line
plotErrors(err_FWD, nominalQ = TRUE) 
plotErrors(err_REV, nominalQ = TRUE)


#### 8) Sample inference ####
# This algorithm will tell us how many unique sequences are in each sample after controlling for sequencing errors. Can also set pool = TRUE to pool across all samples to detect low abundance variants.
dada_FWD <- dada(filt_FWD, err = err_FWD, multithread = TRUE, pool = "pseudo") # took ~ 30 min
dada_REV <- dada(filt_REV, err = err_REV, multithread = TRUE, pool = "pseudo") # took ~ 32 min


#### 9) Merge paired reads ####
# Since we don't have enough overlap to meet the required minimum of 12 bases of overlap (see notes before Step 5), we can lower the minimum overlap with the parameter minOverlap = .  
mergers <- mergePairs(dada_FWD, filt_FWD, dada_REV, filt_REV, verbose = TRUE, minOverlap = 5)
# On average, 92.8% merge success rate across samples

#### 10) Construct a sequence table ####
bact_seq_table <- makeSequenceTable(mergers)
dim(bact_seq_table) 
# (number of samples, number of amplicon sequence variants ASVs) (found 54,916 ASVs across the 57 samples)


#### 11) Remove chimeras ###
bact_seq_table_nochim <- removeBimeraDenovo(bact_seq_table, method = "consensus", multithread = TRUE, verbose = TRUE, minFoldParentOverAbundance = 4)
# found 18,742 chimeras out of 54,916 sequences (34% identified as chimeric)

# We don't know if these chimeras held a lot in terms of abundance. To find out, we can calculate the proportion of sequences retained after chimeras removed.
sum(bact_seq_table_nochim)/sum(bact_seq_table)
# We retained 98% of sequence abundance. Great!


#### 12) Track reads through the pipeline to verify everything worked as expected ####
getN <- function(x) sum(getUniques(x))
(summary_table <- data.frame(row.names = sample.names, 
                             input = out[, 1],
                             filtered = out[, 2], 
                             denoised_FWD = sapply(dada_FWD, getN),
                             denoised_REV = sapply(dada_REV, getN), 
                             merged = sapply(mergers, getN),
                             non_chim = rowSums(bact_seq_table_nochim), 
                             final_perc_reads_retained = round(rowSums(bact_seq_table_nochim)/out[, 1]*100, 1)))
# Looks good! The most amount of reads were removed during filtering, as expected. On average, retained 88.5% of reads.


#### 13) Make summary tables that can be used to assign taxonomy
# give seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(bact_seq_table_nochim)
asv_headers <- vector(dim(bact_seq_table_nochim)[2], mode = "character")

for (i in 1:dim(bact_seq_table_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# make a fasta file of the final ASV seqs (once this is made, upload to rdp.cme.msu.edu/classifier/classifier.jsp)
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/rstudio/Dropbox/SDSU/FG_2022_Bact_ASV.fa")

# make count table
asv_table <- t(bact_seq_table_nochim)
row.names(asv_table) <- sub(">", "", asv_headers)
write.table(asv_table, "FG_2022_Bact_ASV_counts.tsv", sep = "\t", quote = F, col.names = NA)
