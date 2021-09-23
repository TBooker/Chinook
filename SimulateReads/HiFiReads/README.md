## I installed SimLord using the instructions at:
https://bitbucket.org/genomeinformatics/simlord/src/master/

For Julia's tests for building the assembly tutorial, I simulate 30000 long reads. The distribution of read lengths is based on the *C. elegans* model described in the SimLord documentation.

## To activate the conda environment where I can use SimLord...

```
source activate simlord

simlord --read-reference ../../Simulations/Run3/SalmonSim.Stabilising.p1.3.fasta \
        --lognorm-readlength 0.1895 8217 9355 \
        -n 30000 \
        --no-sam \
        --max-passes 2 \
        SalmonSim.Stabilising.p1.3.30k.PacBio

# Zip up that li'l puppy
bgzip  SalmonSim.Stabilising.p1.3.30k.PacBio.fastq

# Then you can minimap2 to align the reads like so:
~/bin/minimap2/minimap2 -ax map-pb ../../ReferenceGenome/SalmonReference.fasta SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz | samtools sort -o SalmonSim.Stabilising.p1.3.30k.PacBio.bam -

samtools index SalmonSim.Stabilising.p1.3.30k.PacBio.bam

## To detach environment...

conda activate base
```

Examining the alignment of these reads, the target coverage of around 40x is hit pretty well.
