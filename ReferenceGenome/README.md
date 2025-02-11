# The simulated reference genome

I downloaded the Chinook salmon reference genome generated by Christensen et al (2018 - PLoS One) from the NCBI database. It was deposited as project *GCA_002872995.1*.

It can be downloaded using:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/872/995/GCF_002872995.1_Otsh_v1.0/GCF_002872995.1_Otsh_v1.0_genomic.fna.gz
```

I extract two 5Mbp regions from this genome to use as a template for the genome I use in SLiM. the first is the firt 5Mbp of Chromosome 30 and the second is from positions 37,500,000 to 42,500,000. The centre of the latter region contains the locus thought to be involved in the red/white colour balanced polymorphism.

The script ```extractFromFasta.py``` makes a reference genome that is appropriate for analysis in SLiM. Note that that script has the name of the genome sequence hard coded in, so if you use something other than Otsh_v1, you'll need to edit the script.

Downstream software requires a file with each chromosome labelled (including the > symbol), that can be easily obtained using:
```
grep ">" SalmonReference.fasta > chrom.list
```
Downstream software requires a file with each chromosome split into FASTA files named for each chromosome:
```
mkdir SalmonReference_split/
gawk -v seq="chr_1" -v RS='>' '$1 == seq {print RS $0}' SalmonReference.fasta > SalmonReference_split/chr_1.fa
gawk -v seq="chr_2" -v RS='>' '$1 == seq {print RS $0}' SalmonReference.fasta > SalmonReference_split/chr_2.fa

```
