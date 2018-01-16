# Calculate the mean and mode of a distribution of protein lengths read in from a file

import sys
from stats import *
from functions import *

BIN_SIZE = 20
MAX_BINS = 1000
binned_data = {}
raw_data = read_list(sys.argv[1])

# Turn into ints
data = []
for r in raw_data:
  data.append(int(r)) 

# Initialize bins
i = 0
while i < MAX_BINS:
  binned_data[(i*BIN_SIZE, (i+1)*BIN_SIZE)] = 0
  i +=1

# Bin the data
for d in data:
  for a, b in binned_data.keys():
    if d > a and d <= b:
      binned_data[(a,b)] += 1
      break

# Output data
max = 0
for a, b in binned_data.keys():
  print "%d\t%d\t%d" % (a, b, binned_data[(a,b)])
