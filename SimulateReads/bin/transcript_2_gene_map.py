## Make transcript to gene mapping files

import sys

transcript_gene = {}

for i in open(sys.argv[1]):
    if len(i) < 10:continue
    gff_attr_str = i.strip().split("\t")[8]
    att_dict = {}
    for a in  gff_attr_str.split(";"):
        att_dict[a.split(" ")[0]] = a.split(" ")[1]
    try:
        transcript_gene[ att_dict["transcript_id"] ] = att_dict["gene"]
    except KeyError:
        continue

for k in transcript_gene.keys():
    print( transcript_gene[k] +" " +k )
