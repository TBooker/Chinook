## A script to generate VCF files for the 

import argparse
import msprime, pyslim, tskit
import numpy as np

## The first thing to do is to read the FASTA file for the reference (there should be only one sequence in this file
## Then, we iterate through the VCF and extract out the haplotypes for an individual, keeping track of the positions
## Then, we do the same with substitutions
## Then, we output a FASTA file with either the two haplotypes or a single string using IUPAC ambiguity codes 


def main():

## Define command line args
	parser = argparse.ArgumentParser(description="This takes result of a SLiM simulation and a refernece genome and outputs the genome of specified individual into FASTA files")

	parser.add_argument("--trees", "-t",
			required = True,
			dest = "trees",
			type = str, 
			help = "The tree sequence output from SLiM.")
			
	parser.add_argument("--output", "-o",
			required = True,
			dest = "vcf",
			type = str, 
			help = "The VCF output from SLiM (or from PySLiM)")
			
	parser.add_argument("--seed",
			required = False,
			dest = "seed",
			type = int, 
			help = "The random seed [default = 666]",
			default = 666)

	args = parser.parse_args()

	sts_tables = pyslim.load( args.trees ).tables

	sts_tables.compute_mutation_times()

	sts = sts_tables.tree_sequence()
#	print( sts.TableCollection )
	
	for v in sts.variants():
		print(v)

	pyslim.SlimTreeSequence( msprime.sim_mutations( sts, rate = 1e-7, random_seed=args.seed ) )
	
	
if __name__ == "__main__":
	main()
