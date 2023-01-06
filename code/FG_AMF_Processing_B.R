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


#### 1) Set working directory ####
path <- "/home/rstudio/Dropbox/SDSU/Faville_Grove_microbial_sequences_2022/AMF_seqs_2022/cutadapt2"
list.files(path)

setwd("/home/rstudio/Dropbox/SDSU")


#### 2) Combine all forward and reverse reads. Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_cut <- sort(list.files(path, pattern ="_R1.fastq.gz", full.names = TRUE))
REV_cut <- sort(list.files(path, pattern ="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
(sample.names <- sapply(strsplit(basename(FWD_cut), "_"), `[`, 2)) #sample ID occurs before the first _ in the fastq file


#### 3) Inspect read quality profiles ####
plotQualityProfile(FWD_cut[54:57])
# These sequences are 300 bp in length
# Quality drops below PHRED score of 30 around 200 bases

plotQualityProfile(REV_cut[54:57])
# Quality drops below PHRED score of 30 around 150 bases


#### 4) Filter and Trim ####
# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Since the primers are already removed, we don't truncate here. Instead, we enforce a minimum sequence length (minLen) to get rid of spurious very low-length sequences. 

FWD_filt <- file.path(path, "filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
REV_filt <- file.path(path, "filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(FWD_filt) <- sample.names
names(REV_filt) <- sample.names

out <- filterAndTrim(FWD_cut, FWD_filt, REV_cut, REV_filt, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Ran into a problem when I checked the first 10 samples; only 12% of reads retained after filtering

