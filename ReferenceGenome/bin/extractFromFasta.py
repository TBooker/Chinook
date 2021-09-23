# A script to extract a region from a FASTA file

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import re
import random

## A quick function to sample nucleotides from A, T, C and G
def repl_fun(match):
	return str(random.choice(list("ATGC")))

## In the Otsh_v1.0 reference genome, NC_037126.1 is name of chromosome 30

for record in SeqIO.parse("ncbi-genomes-2020-09-02/GCF_002872995.1_Otsh_v1.0_genomic.fna", "fasta"):

# We don't care if we are not looking at chr30
	if record.id != "NC_037126.1": continue


# Let's get 10Mbp worth of sequence from the Salmon genome
	start_region = str(record.seq[0:5000000]).upper()

# A 5Mbp region containing the locus thought to be involved in the red/white flesh polymorphism
	pigment_region = str(record.seq[37500000:42500000]).upper()


# Replace any of the IUPAC ambiguity codes or Ns with randomly selected A, C, Ts and Gs
	for l in list('BDEFHIJKLMNOPQRSUVWXYZ'):
		start_region = re.sub(l , repl_fun, start_region)
		pigment_region = re.sub(l , repl_fun, pigment_region)

	region_seq = str(start_region+pigment_region).upper()

	region_seq_record = SeqRecord(
		Seq(region_seq),
		id="SalmonReference",
		description= "")

	start_region_record = SeqRecord(
		Seq(start_region),
		id="chr_1",
		description= "")

	pigment_region_record = SeqRecord(
		Seq(pigment_region),
		id="chr_2",
		description= "")

	SeqIO.write(region_seq_record, "SalmonReferenceForSLiM.fasta", "fasta-2line")
	SeqIO.write([ start_region_record, pigment_region_record ], "SalmonReference.fasta", "fasta")

	break
