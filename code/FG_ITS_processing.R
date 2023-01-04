#### ITS AMPLICON PROCESSING FOR FAVILLE GROVE PROJECT ####
# Soil samples collected from Faville Prairie (remnant) and Tillotson Prairie (restoration) in summer 2022
# Used ITS7 Forward and ITS4 Reverse primers

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


#### 1) Set working directory ####
path <- "/Users/kendallb/Documents/Documents_KK_MacBook_Pro/SDSU_Postdoc/Git/Faville_Grove/data/Microbial_Seqs_2022/FG_ITS_seqs_2022"
list.files(path)

path <- "/home/rstudio/Dropbox/SDSU/Faville_Grove_microbial_sequences_2022/16S_seqs_2022"
list.files(path)

setwd("/home/rstudio/Dropbox/SDSU")


#### 2) Combine all forward and reverse reads. Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_reads <- sort(list.files(path, pattern ="_R1.fastq.gz", full.names = TRUE))
REV_reads <- sort(list.files(path, pattern ="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
(sample.names <- sapply(strsplit(basename(FWD_reads), "_"), `[`, 2)) #sample ID occurs before the first _ in the fastq file


#### 3) Identify that ITS primers are on the FWD and REV reads (ITS7 FWD & ITS4 REV) ####
FWD_primer <- "GTGARTCATCGARTCTTTG"
REV_primer <- "TCCTCCGCTTATTGATATGC"

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


#### 7) Remove primers using cutadapt #### 
# Documentation for cutadapt found at: https://cutadapt.readthedocs.io/en/stable/index.html 
cutadap <- "/Users/kendallb/miniconda3/envs/qiime2-2018.11/bin/cutadapt" # 2.6
cutadapt <- "/Users/kendallb/miniconda3/bin/cutadapt" # 2.6

system2(cutadapt, args = "--version") # Run shell commands from R

# Create output file names for the cutadapt-ed files and create a "cutadapt" folder
path_cut <- file.path(path, "cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)

FWD_cut <- file.path(path_cut, basename(FWD_filtN))
REV_cut <- file.path(path_cut, basename(REV_filtN))

# creat the sequence for the reverse complement of the FWD and REV primers
FWD_revcomp <- dada2:::rc(FWD_primer)
REV_revcomp <- dada2:::rc(REV_primer)

# Trim FWD primer and the reverse-complement of REV primer off of R1 (forward reads)
R1_flags <- paste("-g", FWD_primer, "-a", REV_revcomp) 
# Trim REV and the reverse-complement of FWD primer off of R2 (reverse reads)
R2_flags <- paste("-G", REV_primer, "-A", FWD_revcomp) 

# Run Cutadapt
for(i in seq_along(FWD_reads)) {
  system2(cutadapt, args = c(R1_flags, R2_flags, "-n", 2, 
                             "-o", FWD_cut[i], "-p", REV_cut[i],
                             FWD_filtN[i], REV_filtN[i])) 
}
# -n 2 required to remove FWD and REV primers from reads
# FWD_cut & REV_cut are the cutadapted output files
# FWD_filtN & REV_filtN are the input files

# To verify that the primers were removed, use same code as above earlier, but this time we are looking at the cutadapted reads
rbind(FWD_ForwardReads = sapply(FWD_orients, primerHits, fn = FWD_cut[[1]]), 
      FWD_ReverseReads = sapply(FWD_orients, primerHits, fn = REV_cut[[1]]), 
      REV_ForwardReads = sapply(REV_orients, primerHits, fn = FWD_cut[[1]]), 
      REV_ReverseReads = sapply(REV_orients, primerHits, fn = REV_cut[[1]]))
# Fantastic! The primers are no longer detected in the cutadapted reads


## THE REMAINING STEPS OF THE PIPELINE TAKE PLACE IN THE AWS AMI ##

#### 3) Inspect read quality profiles ####
plotQualityProfile(FWD_cut[1:4])
# These sequences are 250 bp in length
# Quality drops below PHRED score of 30 around 230 bases
# sample Q4 (38th sample) did not sequence properly; only has 3 reads 

plotQualityProfile(REV_cut[1:4])
# Quality drops below PHRED score of 30 around 190 bases
# sample Q4's REV reads don't look good either


#### 4) Filter and Trim ####
# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Since the primers are already removed, we don't truncate here. Instead, we enforce a minimum sequence length (minLen) to get rid of spurious very low-length sequences. 

FWD_filt <- file.path(path, "filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
REV_filt <- file.path(path, "filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(FWD_filt) <- sample.names
names(REV_filt) <- sample.names

out <- filterAndTrim(FWD_cut, FWD_filt, REV_cut, REV_filt, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# On average, 76% of reads were retained after filtering


#### 7) Learn the error rates ####
# Use the filtered reads to learn the sequence error rates and correct for these in the later steps of the pipeline
err_FWD <- learnErrors(FWD_filt, multithread = TRUE) # took less than 2 min
err_REV <- learnErrors(REV_filt, multithread = TRUE) # took less than 3 min

# Visualize the estimated error rates as a sanity check
# red line = expected error rate based on quality score (note that this decreases as the quality increases)
# black line = estimated error rate after convergence of the machine-learning algorithm
# black dots = observed error frequency in our samples
# black dots should align with black line
plotErrors(err_FWD, nominalQ = TRUE) 
plotErrors(err_REV, nominalQ = TRUE)


#### 8) Sample inference ####
# This algorithm will tell us how many unique sequences are in each sample after controlling for sequencing errors. Can also set pool = TRUE to pool across all samples to detect low abundance variants.
dada_FWD <- dada(FWD_filt, err = err_FWD, multithread = TRUE, pool = "pseudo")
dada_REV <- dada(REV_filt, err = err_REV, multithread = TRUE, pool = "pseudo")


#### 9) Merge paired reads ####
# Since we don't have enough overlap to meet the required minimum of 12 bases of overlap (see notes before Step 5), we can lower the minimum overlap with the parameter minOverlap = .  
mergers <- mergePairs(dada_FWD, FWD_filt, dada_REV, REV_filt, verbose = TRUE)
# On average, __ merge success rate across samples


#### 10) Construct a sequence table ####
fungi_seq_table <- makeSequenceTable(mergers)
dim(fungi_seq_table) 
# (number of samples, number of amplicon sequence variants ASVs) (found 7,695 unique sequences across the 57 samples)


#### 11) Remove chimeras ###
fungi_seq_table_nochim <- removeBimeraDenovo(fungi_seq_table, method = "consensus", multithread = TRUE, verbose = TRUE)
# found 31 chimeras out of 7,695 unique sequences (only 0.4% identified as chimeric!)

sum(fungi_seq_table_nochim)/sum(fungi_seq_table)
# retained 99.7% of sequence abundance


#### 12) Track reads through the pipeline to verify everything worked as expected ####
getN <- function(x) sum(getUniques(x))
(summary_table <- data.frame(row.names = sample.names, 
                             input = out[, 1],
                             filtered = out[, 2], 
                             denoised_FWD = sapply(dada_FWD, getN),
                             denoised_REV = sapply(dada_REV, getN), 
                             merged = sapply(mergers, getN),
                             non_chim = rowSums(fungi_seq_table_nochim), 
                             final_perc_reads_retained = round(rowSums(fungi_seq_table_nochim)/out[, 1]*100, 1)))
# Looks ok. The most amount of reads were removed during filtering, as expected. On average, retained 62.4% of reads. 
# Should look back at the mergers step; could be an issue like this one: https://github.com/benjjneb/dada2/issues/213


#### 17) Make summary tables that can be used to assign taxonomy

# giving seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(fungi_seq_table_nochim)
asv_headers <- vector(dim(fungi_seq_table_nochim)[2], mode = "character")

for (i in 1:dim(fungi_seq_table_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}


# make a fasta file of the final ASV seqs (once this is made, upload to rdp.cme.msu.edu/classifier/classifier.jsp)
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/rstudio/Dropbox/SDSU/FG_2022_Fungi_ASV.fa")

# make count table
asv_table <- t(fungi_seq_table_nochim)
row.names(asv_table) <- sub(">", "", asv_headers)
write.table(asv_table, "FG_2022_Fungi_ASV_counts.tsv", sep = "\t", quote = F, col.names = NA)  