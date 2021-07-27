The first thing to do is make a set of SLiM annotations from the CDS in the first 5Mbp and the focal 5Mbp of the Salmon genome

*Make a GFF that is just the sequences annotated as "CDS"*

```
awk '$3 =="CDS"' < GCF_002872995.1_Otsh_v1.0_genomic.gff  > GCF_002872995.1_Otsh_v1.0_genomic.CDS.gff

# Sort the GFF by contig then position within contigs, retaining any header lines
(grep ^"#" GCF_002872995.1_Otsh_v1.0_genomic.CDS.gff; grep -v ^"#" GCF_002872995.1_Otsh_v1.0_genomic.CDS.gff | sort -k1,1 -k4,4n) | bgzip > GCF_002872995.1_Otsh_v1.0_genomic.CDS.sort.gff.gz;

# Generate TABIX index
tabix GCF_002872995.1_Otsh_v1.0_genomic.CDS.sort.gff.gz

# Delete intermediate file
rm GCF_002872995.1_Otsh_v1.0_genomic.CDS.gff.gz

tabix GCF_002872995.1_Otsh_v1.0_genomic.CDS.sort.gff.gz NC_037126.1:1-5000000 > SalmonAnnotation_firstMb.cds.gff

tabix GCF_002872995.1_Otsh_v1.0_genomic.CDS.sort.gff.gz NC_037126.1:37500000-42500000 > SalmonAnnotation_selectedRegion.cds.gff

cat SalmonAnnotation_firstMb.cds.gff SalmonAnnotation_selectedRegion.cds.gff > SalmonAnnotation_forSLiM.gff

# Delete intermediate files:
rm SalmonAnnotation_firstMb.cds.gff
rm SalmonAnnotation_selectedRegion.cds.gff

# merge overlapping elements - we don't care about genome integrity here
~/software/bedtools2/bin/mergeBed -d 2 -i SalmonAnnotation_forSLiM.gff > SalmonAnnotation_forSLiM.merged.bed

 
# Delete intermediate files:
rm SalmonAnnotation_forSLiM.gff

# Convert the output to SLiM readable format 
python BED_to_SLiM.py SalmonAnnotation_forSLiM.merged.bed  > SLiM_annotations_from_GFF.txt

# Delete intermediate files:
rm SalmonAnnotation_forSLiM.merged.bed
```

*Now make a GFF that is all annotated sequences in the relevant regions*

```

# Sort the GFF by contig then position within contigs, retaining any header lines
(grep ^"#" GCF_002872995.1_Otsh_v1.0_genomic.gff; grep -v ^"#" GCF_002872995.1_Otsh_v1.0_genomic.gff | sort -k1,1 -k4,4n) | bgzip > GCF_002872995.1_Otsh_v1.0_genomic.sort.gff.gz;

# Generate TABIX index
tabix GCF_002872995.1_Otsh_v1.0_genomic.sort.gff.gz

tabix GCF_002872995.1_Otsh_v1.0_genomic.sort.gff.gz NC_037126.1:1-5000000 > .gff

python subPositions.py SalmonAnnotation_firstMb.gff 0 > SalmonAnnotation_firstMb.subbed.gff

tabix GCF_002872995.1_Otsh_v1.0_genomic.sort.gff.gz NC_037126.1:37500000-42500000 > SalmonAnnotation_selectedRegion.gff

python subPositions.py SalmonAnnotation_selectedRegion.gff 37500000 > SalmonAnnotation_selectedRegion.subbed.gff

cat SalmonAnnotation_firstMb.subbed.gff SalmonAnnotation_selectedRegion.subbed.gff > SalmonAnnotations_forIGV.gff

# Delete intermediate files:
rm SalmonAnnotation_firstMb.gff
rm SalmonAnnotation_firstMb.subbed.gff
rm SalmonAnnotation_selectedRegion.gff
rm SalmonAnnotation_selectedRegion.subbed.gff

```


The above pipeline will generate a set of annotations for the SLiM script and store it as ```SLiM_annotations_from_GFF.txt```.

I copy and pasted those annotations into the SLiM script I used for the simulations.
