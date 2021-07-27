## Simulate Illumina reads - 
# wgsim is a pretty simplistic read simulator
# the program InSilicoSeq has a but more realism

~/bin/bwa/bwa mem ../../ReferenceGenome/SalmonReference.fasta simulations_sdRAD_1/rad_reads/msp_0.1.fa.gz simulations_sdRAD_1/rad_reads/msp_0.2.fa.gz | samtools sort -o SalmonSim.sdRAD.bam -

samtools index SalmonSim.sdRAD.bam
