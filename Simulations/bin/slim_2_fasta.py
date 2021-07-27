## A script to make make a FASTA file from an individual simulated in SLiM

import argparse
from cyvcf2 import VCF
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

## The first thing to do is to read the FASTA file for the reference (there should be only one sequence in this file
## Then, we iterate through the VCF and extract out the haplotypes for an individual, keeping track of the positions
## Then, we do the same with substitutions
## Then, we output a FASTA file with either the two haplotypes or a single string using IUPAC ambiguity codes 


def main():

## Define command line args
	parser = argparse.ArgumentParser(description="This takes result of a SLiM simulation and a refernece genome and outputs the genome of specified individual into FASTA files")

	parser.add_argument("--ref_fasta", "-r",
			required = True,
			dest = "ref_fasta",
			type = str, 
			help = "The reference genome as a FASTA. This file should only have a single sequence.")
			
	parser.add_argument("--vcf", "-v",
			required = True,
			dest = "vcf",
			type = str, 
			help = "The VCF output from SLiM (or from PySLiM)")

	parser.add_argument("--subs", "-s",
			required = False,
			dest = "subs",
			type = str, 
			help = "The file containing the substitutions")

	parser.add_argument("--ind",
			required = False,
			dest = "ind",
			type = str, 
			help = "The name of the individual in the VCF that you want to extract (the default is the first in the sequence)",
			default = "first")

	parser.add_argument("--diploid",
			required = False,
			dest = "diploid",
			action = "store_true", 
			help = "Use this flag if you want a single reference (this will include IUPAC ambiguity codes)")

	parser.add_argument("--output", "-o",
			required = True,
			dest = "output",
			type = str, 
			help = "The name you want to give to the output file")

	args = parser.parse_args()

## Read the reference genome
## For running SLiM this is a single sequence, I'll split the chromosomes below 

	reference_seqs = list( SeqIO.parse(args.ref_fasta, "fasta") )

	if len(reference_seqs) != 1:
		print("The specified reference genome contains more than one sequence, check your input")
		return
## Get the reference as a single string
	reference = list(reference_seqs[0].seq).copy()

## Get the column index for the individual (the default is to just use the first individual) 
	if args.ind != "first":
		temp_vcf = VCF(args.vcf)
		sample_index = temp_vcf.samples.index(args.ind)
	else:
		sample_index = 0

	substitution_dict = {}
	
## iterate over substitutions and convert them to dict 
	if args.subs:
		for i in open(args.subs):
			if i.startswith("#OUT") or i.startswith("Mutations"): continue
			line = i.strip().split(" ")

# Get the position and nucleotide from the substitutions, then add them to dict
			sub_pos = int(line[3])
			sub_nucl = line[9]

# A quick sanity check
			if sub_nucl not in list("ATGC"): 
				print("What? One of the substitutions is not ATCG")
				return
			substitution_dict[sub_pos -1 ] = [sub_nucl, sub_nucl, "sub"]

## There is a deprecated part of numpy that raises a warning at this stage. I think that it's to do with how CyVCF2 handles the genotype strings

## iterate over variants in the VCF file
	for variant in VCF(args.vcf):
		pos = int(variant.POS)
		
# Add the variants to the substitution dict (this will overwrite any subs that share a position with a variant
		substitution_dict[ pos-1 ] = variant.gt_bases[ sample_index ].split("|") + ["var"]

	iupac_ambiguity = {"AA":"A",
					"CC":"C",
					"GG":"G",
					"TT":"T",
					"AC":"M",
					"AG":"R",
					"AT":"W",
					"CG":"S",
					"CT":"Y",
					"GT":"K",
					"ACG":"V",
					"ACT":"H",
					"AGT":"D",
					"CGT":"B",
					"GAT":"N"}

# Make copies of the reference genome string list
	if args.diploid:
		seq_1 = reference.copy()
	else:
		seq_1 = reference.copy()
		seq_2 = reference.copy()

## Iterate over all substitutions and variants and substitute the appropriate strings into the ref_seqs
	for p in substitution_dict.keys():

		bases = substitution_dict[p][:2]
## If you want a diploid sequence, 
		if args.diploid:
			seq_1[p] = iupac_ambiguity["".join(sorted( bases)) ]

		else:
			seq_1[p] = bases[0]
		
			seq_2[p] = bases[1]
		
# Put the sequence(s) into a iterable list
	if args.diploid:
		raw_seqs = ["".join( seq_1 )]
	else:
		raw_seqs = ["".join( seq_1 ), "".join( seq_2 )]

	output_seqs = []
	
	for n in range(len(raw_seqs)):

		hap = str( n + 1 ) ## make a 1-based label to add to the file

		chrom_1 = raw_seqs[n][:5000000] # Each of the simulated chromosomes were 5Mbp long
		chrom_2 = raw_seqs[n][5000000:]

		chrom_1_record = SeqRecord(
			Seq(chrom_1),
			id="chr_1_"+hap,
			description= "")

		chrom_2_record = SeqRecord(
			Seq(chrom_2),
			id="chr_2_"+hap,
			description= "")

		output_seqs.append( chrom_1_record )
		output_seqs.append( chrom_2_record )

	SeqIO.write(output_seqs, args.output, "fasta")

if __name__ == "__main__":
	main()
