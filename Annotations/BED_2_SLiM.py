# Quick script to get the CDS annotations from Salmon into a format that can be pasted into the SLiM config file

import sys, gzip
## I merged the GFF annotations

last_end = 0
length = 0
for i in open(sys.argv[1]):
	x = i.strip().split()
	start = int(x[1]) 
	end = int(x[2])
	if start > 37.5e6:
		start = (start - 37500000)  + 5000000
		end = (end - 37500000) + 5000000
#	print( "initializeGenomicElement(g1, " + str(last_end) +", " + str(start - 1) + ");" )  # neutral spacers
	print( "initializeGenomicElement(g2, " + str(start) +", " + str(end) + ");" ) # CDS
	length += end- start

	last_end = end+1
#print( "initializeGenomicElement(g1, " + str(last_end) +", " + str(10000000-1) + ");" )



