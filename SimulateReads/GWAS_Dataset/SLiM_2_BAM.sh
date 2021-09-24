## Simulate Illumina reads -
# wgsim is a pretty simplistic read simulator
# the program InSilicoSeq has a but more realism

# 800,000 corresponds to around 10x coverage - for GWAS data
r=800000

#for pop in 1 3 5 7 9
for pop in 1 3
  do
#    for ind in $(seq 1 5)
    for ind in $(seq 1 2)
      do
        python ../../Simulations/bin/slim_2_fasta.py --ref_fasta ../../ReferenceGenome/SalmonReferenceForSLiM.fasta --vcf ../../Simulations/Run3/SalmonSim.Stabilising.p$pop.3.vcf --subs ../../Simulations/Run3/SalmonSim.Stabilising.Substitutions.txt  --ind i$ind --output SalmonSim.Stabilising.p$pop.3.i$ind.fa
        iss generate -g SalmonSim.Stabilising.p$pop.3.i$ind.fa --model hiseq --output SalmonSim.Stabilising.p$pop.3.i$ind.$r --cpus 6 --n_reads $r

        bgzip SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R1.fastq
        bgzip SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R2.fastq

        ~/bin/bwa/bwa mem ../../ReferenceGenome/SalmonReference.fasta SalmonSim.Stabilising.p${pop}.3.i${ind}.${r}_R1.fastq.gz SalmonSim.Stabilising.p${pop}.3.i${ind}.${r}_R2.fastq.gz -t4 | samtools sort -o Chinook.p${pop}.i${ind}.bam -

        samtools index Chinook.p${pop}.i${ind}.bam

# Clean up intermediates
        rm SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R1.fastq.gz
        rm SalmonSim.Stabilising.p$pop.3.i$ind.${r}_R2.fastq.gz
        rm SalmonSim.Stabilising.p$pop.3.i$ind.fa

      done
    done
