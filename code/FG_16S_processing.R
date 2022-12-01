#### 16S AMPLICON PROCESSING FOR FAVILLE GROVE PROJECT ####
# Soil samples collected from Faville Prairie (remnant) and Tillotson Prairie (restoration) in summer 2022
# Used the Earth Microbiome Project primers: 515F/806R

(.packages())

library(dada2)
packageVersion("dada2")

install.packages("Biostrings")

#### 1) Set working directory ####
path <- "/Users/kendallb/Documents/Documents_KK_MacBook_Pro/SDSU_Postdoc/Faville_Grove_prj/Faville_Grove_Seqs_2022/Faville_Grove_16S_seqs_2022"
list.files(path)

path <- "/Users/kendallb/Documents/Documents_KK_MacBook_Pro/SDSU_Postdoc/Git/Faville_Grove/data/Microbial_Seqs_2022/FG_16S_seqs_2022"
list.files(path)

#### 2) Combine all forward and reverse reads. Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_reads <- sort(list.files(path, pattern ="_R1.fastq.gz", full.names = TRUE))
REV_reads <- sort(list.files(path, pattern ="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
(sample.names <- sapply(strsplit(basename(FWD_reads), "_"), `[`, 1)) #sample ID occurs before the first _ in the fastq file


#### 3) Inspect read quality profiles ####
plotQualityProfile(FWD_reads[3:6])
# Quality doesn't drop below q30; no need to truncate heavily

plotQualityProfile(REV_reads[3:6])
# Quality doesn't drop below q30; no need to truncate heavily

# F primer: GTGYCAGCMGCCGCGGTAA (spans base positions 515-533)
# R primer: GGACTACNVGGGTWTCTAA (spans base positions 806-824)
# Target amplicon sequence is 274 bps in length (806-533 + 1)
# We have 300 bps of length with the F & R reads combined (150 bps each), so 300 bps - 274 (target amplicon length) = 26 bps of overlap
# dada2's mergePairs command needs at least 12 bases to be overlapping, but let's say 20 to be conservative. So, 26-20 = 6; 6 divided by 2 = 3 bps could be trimmed from each F & R read
# Ben Callahan recommends adding 15 bases to the length of the target amplicon; truncation lengths of F & R reads should sum to this number, as a minimum: 274 + 15 = 289; If I truncate to 148 bases on both F & R reads, that sums to 296 

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


#### 5) Filter and Trim ####
# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note that the primers are still are on the forward & reverse reads, so we remove them with the trimLeft parameter(length of FWD primer, length of REV primer). We used the 341 forward primer: CCTACGGGNGGCWGCAG (17 bp in length) & 785 reverse primer: GACTACHVGGGTATCTAATCC (21 bps in length).
filt_FWD <- file.path(path, "filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
filt_REV <- file.path(path, "filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(filt_FWD) <- sample.names
names(filt_REV) <- sample.names
out <- filterAndTrim(FWD_reads, filt_FWD, REV_reads, filt_REV, trimLeft = c(19, 19), truncLen = c(148, 148),
                     maxN = 0, maxEE = c(2,5), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE, verbose = TRUE)