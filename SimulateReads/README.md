# Simulate Reads

This is the core of the whole thing, using the simulated genomes to generate sequence data.

All the software used here is really easy to use off-the-shelf. I've noted a couple of quirks I came across along the way, but mostly things were straightforward.

Generally, we are working with a 10Mbp genome, so we can use the Lander-Waterman equation to figure out how many reads we'll need to get a target coverage. For example, if we wanted a mean coverage of 10x using paired-end Illumina HiSeq reads, which have an average length of 125bp per end, we would need:

<img src="https://latex.codecogs.com/gif.latex?Number of Reads = \frac{(Target Coverage \times Genome Size)}{Sequence Length} " />

Plugging the numbers in, we'd need around 800000 reads for an average of 10x coverage in that case.


# Illumina Reads

We use the simulator InSilicoSeq for generating Illumina-like reads.

The only hitch with InSilicoSeq was that there seemed to be a conflict with the version of BioPython I had been using and I ended up having to downgrade to an earlier version (see issue #208 on the InSilicoSeq Github page).

The script ```Illumina_ISS/SLiM_2_BAM.sh``` takes a single argument as input, the fasta file of an individual's genome. That is then used to generate the reads we require.


# Long reads

We use the simulator SimLord

The only hitch with SimLord is that the authors recommend using a particular conda environment, which you have to attach before using and detach when you're done. was that there seemed to be a conflict with the version of BioPython I had been using and I ended up having to downgrade to an earlier version (see issue #208 on the InSilicoSeq Github page).

The script ```Illumina_ISS/SLiM_2_BAM.sh``` takes a single argument as input, the fasta file of an individual's genome. That is then used to generate the reads we require.
