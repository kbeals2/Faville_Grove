#### 16S AMPLICON PROCESSING FOR FAVILLE GROVE PROJECT ####
# Soil samples collected from Faville Prairie (remnant) and Tillotson Prairie (restoration) in summer 2022
# Used the Earth Microbiome Project primers: 515F/806R

(.packages())

install.packages("BiocManager")
library(BiocManager)

BiocManager::install("Biostrings")
library(Biostrings)

BiocManager::install("ShortRead")
library(ShortRead)

library(dada2)
packageVersion("dada2")

#### 1) Set working directory ####
path <- "/Users/kendallb/Documents/Documents_KK_MacBook_Pro/SDSU_Postdoc/Git/Faville_Grove/data/Microbial_Seqs_2022/FG_16S_seqs_2022"
list.files(path)

#### 2) Combine all forward and reverse reads. Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_reads <- sort(list.files(path, pattern ="_R1.fastq.gz", full.names = TRUE))
REV_reads <- sort(list.files(path, pattern ="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
(sample.names <- sapply(strsplit(basename(FWD_reads), "_"), `[`, 1)) #sample ID occurs before the first _ in the fastq file


#### 3) Inspect read quality profiles ####
plotQualityProfile(FWD_reads[1:4])
# Quality doesn't drop below q30; no need to truncate heavily

plotQualityProfile(REV_reads[1:4])
# Quality doesn't drop below q30; no need to truncate heavily


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

# The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. We are going to “pre-filter” the sequences just to remove those with Ns.
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
# Target amplicon sequence is 274 bps in length (806-533 + 1)
# We have 300 bps of length with the F & R reads combined (150 bps each), so 300 bps - 274 (target amplicon length) = 26 bps of overlap
# dada2's mergePairs command needs at least 12 bases to be overlapping, but let's say 20 to be conservative. So, 26-20 = 6; 6 divided by 2 = 3 bps could be trimmed from each F & R read
# Ben Callahan recommends that the sum of the truncated FWD & REV reads should be 20 bases longer than the length of the sequenced amplicon. If the reads are truncated to 149 bases each, that sums to 298 bases. 274 bases (length of the amplicon) + 20 = 294.


#### 5) Filter and Trim ####
# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note that the primers are still on the FWD & REV reads, so we remove them with the trimLeft parameter(length of FWD primer = 19 bases, length of REV primer = 19 bases).
filt_FWD <- file.path(path, "filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
filt_REV <- file.path(path, "filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(filt_FWD) <- sample.names
names(filt_REV) <- sample.names
out <- filterAndTrim(FWD_reads, filt_FWD, REV_reads, filt_REV, trimLeft = c(19, 19), truncLen = c(149, 149),
                     maxN = 0, maxEE = c(2,5), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE, verbose = TRUE)


#### 6) Inspect read quality profiles after trimming ####
plotQualityProfile(filt_FWD[1:4])
plotQualityProfile(filt_REV[1:4])


#### 7) Learn the error rates ####
# Use the filtered reads to learn the sequence error rates and correct for these in the later steps of the pipeline
err_FWD <- learnErrors(filt_FWD, multithread = TRUE) # took ~ 11 min
err_REV <- learnErrors(filt_REV, multithread = TRUE) # took ~ 16 min

# Visualize the estimated error rates as a sanity check
# red line = expected error rate based on quality score (note that this decreases as the quality increases)
# black line = estimated error rate after convergence of the machine-learning algorithm
# black dots = observed error frequency in our samples
# black dots should align with black line
plotErrors(err_FWD, nominalQ = TRUE) 
plotErrors(err_REV, nominalQ = TRUE)


#### 8) Sample inference ####
# his algorithm will tell us how many unique sequences are in each sample after controlling for sequencing errors. Can also set pool = TRUE to pool across all samples to detect low abundance variants
dadaFs <- dada(filt_FWD, err = err_FWD, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filt_REV, err = err_REV, multithread = TRUE, pool = "pseudo")


#### 9) Merge paired reads ####
mergers <- mergePairs(dadaFs, filt_FWD, dadaRs, filt_REV, verbose = TRUE)
