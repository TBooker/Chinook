## Simulate Illumina reads -
# wgsim is a pretty simplistic read simulator
# the program InSilicoSeq has a but more realism

# 800,000 corresponds to around 10x coverage - for IGV demo
# 1,600,000 corresponds to around 20x coverage - for IGV demo
# 6,400,000 corresponds to around 80x coverage - for genome assembly


for r in 800000 1600000 6400000
do
    iss generate -g ../../Simulations/Run3/SalmonSim.Stabilising.p1.3.fasta --model hiseq --output SalmonSim.Stabilising.p1.1.$r --cpus 6 --n_reads $r

    bgzip SalmonSim.Stabilising.p1.1.${r}_R1.fastq
    bgzip SalmonSim.Stabilising.p1.1.${r}_R2.fastq

    ~/bin/bwa/bwa mem ../../ReferenceGenome/SalmonReference.fasta SalmonSim.Stabilising.p1.1.${r}_R1.fastq.gz SalmonSim.Stabilising.p1.1.${r}_R2.fastq.gz | samtools sort -o SalmonSim.Stabilising.p1.1.$r.bam -

    samtools index SalmonSim.Stabilising.p1.1.$r.bam

#    ~/bin/bcftools/bcftools mpileup -f ../../ReferenceGenome/SalmonReference.fasta SalmonSim.Stabilising.p1.1.$r.bam | ~/bin/bcftools/bcftools call -mv > SalmonSim.Stabilising.p1.1.$r.raw.vcf

#    ~/bin/bcftools/bcftools filter -s LowQual -e '%QUAL<30 || DP>100' SalmonSim.Stabilising.p1.1.$r.raw.vcf  > SalmonSim.Stabilising.p1.1.$r.raw.vcf
done
