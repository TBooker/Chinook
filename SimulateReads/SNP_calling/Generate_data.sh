## Simulate Illumina reads for the GWAS tutorial

# Make locations for file outputs
mkdir gvcf
mkdir db
mkdir vcf
mkdir log
mkdir abundance

# Set locations of programs
picard=~/bin/picard.jar
gatk=~/bin/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar


# 400,000 corresponds to around 10x coverage - for GWAS data
r=400000

# Make a directory to write the BAM files into
mkdir bam
for r in 400000 800000 1600000
  do
    for pop in  1 10
    do
      for ind in 1 2
        do
          python ../../Simulations/bin/slim_2_fasta.py --ref_fasta ../../ReferenceGenome/SalmonReferenceForSLiM.fasta --vcf ../../Simulations/Run3/SalmonSim.Stabilising.p$pop.3.vcf --subs ../../Simulations/Run3/SalmonSim.Stabilising.Substitutions.txt  --ind i$ind --output SalmonSim.Stabilising.p$pop.3.i$ind.fa
          iss generate -g SalmonSim.Stabilising.p$pop.3.i$ind.fa --model hiseq --output Chinook.p$pop.3.i$ind.$r --cpus 6 --n_reads $r

          bgzip Chinook.p$pop.3.i$ind.${r}_R1.fastq
          bgzip Chinook.p$pop.3.i$ind.${r}_R2.fastq

# Align reads to the actual genome seuqence used...
          ~/bin/bwa/bwa mem ../../ReferenceGenome/SalmonReference.fasta Chinook.p${pop}.3.i${ind}.${r}_R1.fastq.gz Chinook.p${pop}.3.i${ind}.${r}_R2.fastq.gz -t4 | samtools sort -o bam/Chinook.p${pop}.i${ind}.r${r}.bam -


          java -jar $picard  AddOrReplaceReadGroups \
                I=bam/Chinook.p${pop}.i${ind}.r${r}.bam \
                O= bam/Chinook.p${pop}.i${ind}.r${r}.rg.bam \
                RGID=1 \
                RGLB=lib1 \
                RGPL=illumina \
                RGPU=unit1 \
                RGSM=20

            samtools index bam/Chinook.p${pop}.i${ind}.r${r}.rg.bam


# Clean up intermediates
        rm bam/Chinook.p${pop}.i${ind}.bam
        rm SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R1.fastq.gz
        rm SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R2.fastq.gz
        rm SalmonSim.Stabilising.p$pop.3.i$ind.fa

        mv *_abundance.txt abundance/

      done

done
