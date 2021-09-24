#!/usr/bin/env Rscript

## This script generates the data for the RNA-seq tutorial for BIOL525D
## It uses the reference genome from Chinook, but I don't actually use
## the population data to motivate the expression differences

args <- commandArgs(trailingOnly = TRUE)

if (length(args) <3){
  stop("Usage: Rscript --vanilla reference_genome_loc salmon_gtf_file outdir")
}

# I lifted this gffRead function  from the ballgown package...
# It reads in a GFF file and labels the columns appropriately

gffRead = function (gffFile, nrows = -1, verbose=FALSE) 
{
  if(verbose){
    cat("Reading ", gffFile, ": ", sep = "")    
  }
  gff = read.table(gffFile, sep = "\t", as.is = TRUE, quote = "", 
                   header = FALSE, comment.char = "#", nrows = nrows, 
                   colClasses = c("character", "character", "character", "integer", 
                                  "integer", "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", 
                    "end", "score", "strand", "frame", "attributes")
  if(verbose){
    cat("found", nrow(gff), "rows with classes:", paste(sapply(gff, 
                                                               class), collapse = ", "), "\n")
  }
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

## Here I'll use Polyester to generate RNA seq reads from the Chinook simulations

## Install...
# source("http://bioconductor.org/biocLite.R")
# biocLite("polyester")
## Commented out cause I already did it.
## I developed this using R 3.4, for R >4, the biocLite installation function is different 


## We have 414 transcripts in the Chinook genome
## Let's specify the numbers of reads to get from each...
## Add some noise to the read counts
# This is less preferable than using a model based approach where gene length dictates the num of reads
tx_lengths = read.csv("~/work/Chinook/SimulateReads/RNA/dev/transcript_lengths.txt", header = F, sep = " ")
#Let's normalise the transcript lengths and use that to inform read numbers.
# the following gives reasonable numbers 
read_counts = as.integer(500 + 300*(tx_lengths$V2 - mean(tx_lengths$V2))/sd(tx_lengths$V2))
#read_counts = as.integer(rnorm(414, 1, 0.1)* 300)

# Here's the gene expression differences between treatments - give it random noise using the Poisson dist
fold_change_mat = matrix(c(rpois(414, 1)+1, rpois(414, 1)+1), nrow = 414)

# Make 10 of the genes be up regulated in each population...
# - draw the up and down regulation from a Poisson Distribution
fold_change_mat[sample(414,10),1] = rpois(10,10)
fold_change_mat[sample(414,10),2] = rpois(10,10)


library(polyester)
library(Biostrings)

# Specify the location of the Chinook reference genome
# polyester requires a dir containing each chromosome in a separate FASTA
#reference_genome_loc = "~/work/Chinook/ReferenceGenome/SalmonReference_split/"
reference_genome_loc = args[1]

# Read in the Salmon GFF - polyester uses this to build the reads
#salmon_gtf_file = "~/work/Chinook/Annotations/SalmonAnnotations_forPolyester.gff"
salmon_gtf_file = args[2]

gtf = gffRead(salmon_gtf_file)

simulate_experiment(gtf = gtf, 
                    seqpath = reference_genome_loc,
        # The number of biolofical replicates to sim for each column of the expression matrix 
                    num_reps=c(3, 3),
        # How do you want to specify the number of reads for each transcript?            
#                    meanmodel = TRUE, 
                    reads_per_transcript = read_counts,
        # Specify the expresssion matrix
                    fold_changes = fold_change_mat, 
        # Specify the output location
                    outdir=args[3], 
        # How are attributes separated in the GFF?
                    attrsep= ";",
        # What error model would you like to use for the reads?
                    error_model = "illumina5") 





