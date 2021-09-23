# Make a little file containing the lengths of transcripts

import sys

tx_dict = {}
summer = 0
for i in open(sys.argv[1]):
  x =i.strip().split(" ")
  tx_dict[ x[0].split("na-")[1] ] = int(x[1])
  summer+=int(x[1])
#print( len(tx_dict.keys()))
mean_length = int(summer/len(tx_dict.keys()))

for j in open(sys.argv[2]):
    x =j.strip()
    try:
        print(x, tx_dict[x])
    except KeyError:
        print(x, str(mean_length))
