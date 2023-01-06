#### ITS AMPLICON PROCESSING FOR FAVILLE GROVE PROJECT ####
# Soil samples collected from Faville Prairie (remnant) and Tillotson Prairie (restoration) in summer 2022
# Used ITS7 Forward and ITS4 Reverse primers

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

#### 1) Set working directory ####
path <- "/Users/kendallb/Documents/Documents_KK_MacBook_Pro/SDSU_Postdoc/Git/Faville_Grove/data/Microbial_Seqs_2022/FG_ITS_seqs_2022"
list.files(path)


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


#### 4) Remove primers using cutadapt #### 
# Documentation for cutadapt found at: https://cutadapt.readthedocs.io/en/stable/index.html 
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
