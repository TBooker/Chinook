## Simulate Illumina reads - 
# wgsim is a pretty simplistic read simulator
# the program InSilicoSeq has a but more realism

for r in 500000 2000000
do
    iss generate -g ../../Simulations/Run3/SalmonSim.Stabilising.p1.3.fasta --model hiseq --output SalmonSim.Stabilising.p1.1.$r --cpus 6 --n_reads $r

    ~/bin/bwa/bwa mem ../../ReferenceGenome/SalmonReference.fasta SalmonSim.Stabilising.p1.1.${r}_R1.fastq SalmonSim.Stabilising.p1.1.${r}_R2.fastq | samtools sort -o SalmonSim.Stabilising.p1.1.$r.bam -

    samtools index SalmonSim.Stabilising.p1.1.$r.bam

    ~/bin/bcftools/bcftools mpileup -f ../../ReferenceGenome/SalmonReference.fasta SalmonSim.Stabilising.p1.1.$r.bam | ~/bin/bcftools/bcftools call -mv > SalmonSim.Stabilising.p1.1.$r.raw.vcf

    ~/bin/bcftools/bcftools filter -s LowQual -e '%QUAL<30 || DP>100' SalmonSim.Stabilising.p1.1.$r.raw.vcf  > SalmonSim.Stabilising.p1.1.$r.raw.vcf
done
