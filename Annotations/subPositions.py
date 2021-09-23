import sys

if len(sys.argv) < 3:
	print("USAGE:\n\tpython subPositions.py InputGFF numberOfBases")
	sys.exit()

for i in open(sys.argv[1]):
	line = i.strip().split("\t")
	if line[2] == "region": continue

	if int(line[3]) <= 5000000 or int(line[4]) <= 5000000:
		line[0]= "chr_1"
	elif int(line[3]) > 5000000:
		line[0]= "chr_2"
#	if not line[0].startswith("c"): print("!!!")
	if int(line[3]) < int(sys.argv[2]): continue
	start = str( int(line[3])- int(sys.argv[2]) )
	stop = str( int(line[4]) - int(sys.argv[2]) )
	line[3] = start
	line[4] = stop
	print( "\t".join(line) )
