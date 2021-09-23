# Simulate an RNA seq experiment to generate data for the tutorial

Here I use the Salmon reference genome to simulate an RNA-seq experiment. In this hypothetical experiment, let's pretend that we are comparing fish from up and downstream of the river we simulated. The up stream populations may experience cooler temperatures due to snowmelt. Along with all the local adaptation that we can see between these populations, there may be gene expression differences between the different locations.

The script ```generateRNAseqReads.R``` uses the program "polyester" (Frazee et al 2015) to simulate sequencing reads for an RNA seq experiment. We simulate an RNA seq experiment consisting of 3 individuals from upstream (cold) and 3 individuals form downstream (hot).

Here's how I invoke it:


After a preliminary run, I had a reference containing the transcript sequences for each gene as a FASTA. I used that to calculate transcript length and to specify the number of reads in polyester.

```
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next} {print val,length($0)}' SalmonReference_ref.transcripts.fa  > TranscriptLengths.txt
```

Polyester does not simulate varying read quality, but that's really not of central importance to this part of the tutorial so let's ignore that and just add phoney, high read qualities to the FASTA files that polyester generates to make FASTQ files. I use a little script called ```fasta_2_fastq.py``` for this:

```

Rscript --vanilla  generateRNAseqReads.R ../../ReferenceGenome/SalmonReference_split/ ../../Annotations/SalmonAnnotations_forPolyester.gff Warm_Cold_Experiment

parallel "python ../bin/fasta_2_fastq.py Warm_Cold_Experiment/sample_{1}_{2}.fasta > Warm_Cold_Experiment/sample_{1}_{2}.fq" ::: 01 02 03 04 05 06 ::: 1 2
gzip Warm_Cold_Experiment/*fq
rm Warm_Cold_Experiment/*fasta

parallel "mv Warm_Cold_Experiment/sample_{1}_{2}.fq.gz  Warm_Cold_Experiment/warm_sample_{1}_{2}.fq.gz" ::: 01 02 03 ::: 1 2

parallel "mv Warm_Cold_Experiment/sample_{1}_{2}.fq.gz  Warm_Cold_Experiment/cold_sample_{1}_{2}.fq.gz" ::: 04 05 06 ::: 1 2

```
The file ```Warm_Cold_Experiment/sim_tx_info.txt``` contains the actual fold change in expression that was simulated, so provides the ground truth against which to the compare the output of the analysis.



## Below are the steps to analyse gene expression using RSEM - following the BIOL525D tutorial

I made a few changes to the tutorial that Kay developed. Specifically, we do not start with a set of transcripts. We start with a GFF file and a reference genome. From those we make the transcript reference against which we align the RNA-seq reads.

Here're the first steps:

```
### Make the RSEM reference

mkdir TutorialDemo
cd TutorialDemo

rsem-prepare-reference --gff3 ../../../Annotations/SalmonAnnotations_forIGV.gff ../../../ReferenceGenome/SalmonReference.fasta SalmonReference_ref

# Make the BowTie Reference

bowtie2-build -f SalmonReference_ref.transcripts.fa SalmonReference_ref
```


### Use RSEM to map the reads using BowTie and then quantify expression:

```
# For a single sample
rsem-calculate-expression --bowtie2 --paired-end ../Warm_Cold_Experiment/warm_sample_01_1.fq.gz ../Warm_Cold_Experiment/warm_sample_01_2.fq.gz SalmonReference_ref warm_sample_01

## Takes about 1-2 mins per sample

## Doing this for the other samples in parallel:

parallel "rsem-calculate-expression --bowtie2 --paired-end ../Warm_Cold_Experiment/{1}_sample_{2}_1.fq.gz ../Warm_Cold_Experiment/{1}_sample_{2}_2.fq.gz SalmonReference_ref {1}_sample_{2}" ::: warm ::: 02 03
    ## missing 01 as we already did that above

parallel "rsem-calculate-expression --bowtie2 --paired-end ../Warm_Cold_Experiment/{1}_sample_{2}_1.fq.gz ../Warm_Cold_Experiment/{1}_sample_{2}_2.fq.gz SalmonReference_ref {1}_sample_{2}" ::: cold ::: 04 05 06

cd ../

```
