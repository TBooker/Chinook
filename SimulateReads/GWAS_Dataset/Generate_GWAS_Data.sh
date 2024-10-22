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

for pop in $(seq 1 10)
  do
    for ind in $(seq 0 9)
#    for ind in $(seq 1 2)
      do
        python ../../Simulations/bin/slim_2_fasta.py --ref_fasta ../../ReferenceGenome/SalmonReferenceForSLiM.fasta --vcf ../../Simulations/Run3/SalmonSim.Stabilising.p$pop.3.vcf --subs ../../Simulations/Run3/SalmonSim.Stabilising.Substitutions.txt  --ind i$ind --output SalmonSim.Stabilising.p$pop.3.i$ind.fa
        iss generate -g SalmonSim.Stabilising.p$pop.3.i$ind.fa --model hiseq --output SalmonSim.Stabilising.p$pop.3.i$ind.$r --cpus 6 --n_reads $r

        bgzip SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R1.fastq
        bgzip SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R2.fastq

# Align reads to the actual genome seuqence used...
        ~/bin/bwa/bwa mem ../../ReferenceGenome/SalmonReference.fasta SalmonSim.Stabilising.p${pop}.3.i${ind}.${r}_R1.fastq.gz SalmonSim.Stabilising.p${pop}.3.i${ind}.${r}_R2.fastq.gz -t4 | samtools sort -o bam/Chinook.p${pop}.i${ind}.bam -


        java -jar $picard  AddOrReplaceReadGroups \
              I=bam/Chinook.p${pop}.i${ind}.bam \
              O=bam/Chinook.p${pop}.i${ind}.rg.bam \
              RGID=1 \
              RGLB=lib1 \
              RGPL=illumina \
              RGPU=unit1 \
              RGSM=20

        samtools index bam/Chinook.p${pop}.i${ind}.rg.bam


# Clean up intermediates
        rm bam/Chinook.p${pop}.i${ind}.bam
        rm SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R1.fastq.gz
        rm SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R2.fastq.gz
        rm SalmonSim.Stabilising.p$pop.3.i$ind.fa

        mv *_abundance.txt abundance/

      done

done

## Now here's the code from Julia's SNP calling tutorial - I use this to make the files consistant with what the students will be working with...
ls bam/ | grep .rg.bam$ | sed s/.rg.bam//g > samplelist.txt

cat samplelist.txt


while read name; do
  java -jar $picard MarkDuplicates \
  I=bam/$name.rg.bam O=bam/$name.dedup.bam \
  M=log/$name.duplicateinfo.txt
  samtools index bam/$name.dedup.bam
  rm bam/$name.rg.bam
  rm bam/$name.rg.bam.bai
done < samplelist.txt

java -jar $picard CreateSequenceDictionary R= ../../ReferenceGenome/SalmonReference.fasta O= ../../ReferenceGenome/SalmonReference.dict

samtools faidx ../../ReferenceGenome/SalmonReference.fasta

#gatk-package-4.2.2.0-local.jar
for name in `cat samplelist.txt`
  do
    java -Xmx15g -jar $gatk HaplotypeCaller \
    -R ../../ReferenceGenome/SalmonReference.fasta \
    -I bam/$name.dedup.bam \
    --native-pair-hmm-threads 6 \
    -ERC GVCF \
    -O gvcf/$name.dedup.g.vcf
  done


# Make sample map for GATK
for i in `ls gvcf/*g.vcf | sed 's/.dedup.g.vcf//g' | sed 's/gvcf\///g'`
do
  echo -e "$i\tgvcf/$i.dedup.g.vcf"
done > Chinook_GWAS.sample_map

# Build the DB for the samples and then call a single VCF for all individuals - one VCF per chromosome
for i in chr_1 chr_2
  do
    java -Xmx10g -Xms10g -jar $gatk \
       GenomicsDBImport \
       --genomicsdb-workspace-path db/SalmonReference_${chr} \
       --batch-size 50 \
       -L chr_1 \
       --sample-name-map Chinook_GWAS.sample_map \
       --reader-threads 3

    java -Xmx10g -jar $gatk GenotypeGVCFs \
      -R ../../ReferenceGenome/SalmonReference.fasta \
      -V gendb://db/SalmonReference_${chr} \
      -O vcf/Chinook_GWAS.${chr}.vcf.gz
    done

# Concatenate the VCFs for each chromosome into a single file...
bcftools concat \
  vcf/Chinook_GWAS.chr_1.vcf.gz \
  vcf/Chinook_GWAS.chr_2.vcf.gz \
  -O z > vcf/Chinook_GWAS.vcf.gz
