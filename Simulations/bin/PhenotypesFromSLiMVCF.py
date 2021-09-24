
import argparse
from cyvcf2 import VCF
import numpy as np

def main():
## Define command line args
    parser = argparse.ArgumentParser(description="This takes result of a SLiM simulation and a refernece genome and outputs the genome of specified individual into FASTA files")

    parser.add_argument("--vcf", "-v",
            required = True,
            dest = "vcf",
            type = str,
            help = "The vcf file you want to use.")

    parser.add_argument("--output", "-o",
            required = True,
            dest = "output",
            type = str,
            help = "The name of the output file")

    args = parser.parse_args()

    phenotypic_contributions = []
## iterate over variants in the VCF file
    for variant in VCF(args.vcf):
        if variant.INFO["S"] == 0:continue
        pos = int(variant.POS)

        if type(variant.INFO["S"] ) == float:
            selCoeffs = np.array([0., variant.INFO["S"]])
        elif type(variant.INFO["S"] ) == tuple:
            selCoeffs = np.array([0] + list(variant.INFO["S"]))
        if selCoeffs.sum() == 0: continue
        phenotypic_contribution_from_variant = []
        for g in  variant.genotypes:
            phenotypic_contribution_from_variant.append((sum( [selCoeffs[t] for t in g[:-1]] ) ))
        phenotypic_contributions.append(phenotypic_contribution_from_variant)

    phenotype_array =  np.array( phenotypic_contributions )
    phenotypes =  phenotype_array.sum(axis = 0)

    pop = args.vcf.split("/")[-1].split(".")[2]

    output_file = open( args.output, "w")
    output_file.write("Individual,Phenotype\n")

    for ind,phen  in zip(range(len(phenotypes)), phenotypes):
        output_file.write("Chinook."+pop+".i"+str(ind) + ","+ str(phen)+"\n")
#    Chinook.p8.i8

main()
