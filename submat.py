# Read in two lists
#  1. ordered list of two-letter organism codes e.g. ['ab', 'cc']
#  2. ordered list of Pfam domains e.g. ['7tm_1', 'AAA']


# Read in a matrix looking like
#  	ab 	cc 	cn 	cp
# 7tm_1	0 	0 	1 	0
# AAA	38 	43 	31 	47
# ABC_tran	36 	46 	34

# Output a submatrix...

import sys
from functions import read_matrix, read_list

TRUE = 1
FALSE = 0

OMIT_ZERO_ROWS = FALSE  # Whether to leave out rows with only zeros

try:
  dat_path = sys.argv[1]
  rows_path = sys.argv[2]
  cols_path = sys.argv[3]
except:
  print "ARGS: <matrix file> <list of rows> <list of columns>"
  sys.exit()

sys.stderr.write("Reading in data matrix %s\n" % dat_path)
D = read_matrix(dat_path)

sys.stderr.write("Reading in rows%s\n" % rows_path)
rows = read_list(rows_path)
sys.stderr.write("%d rows\n" % len(rows))

sys.stderr.write("Reading in columns %s\n" % cols_path)
cols = read_list(cols_path)
sys.stderr.write("%d columns\n" % len(cols))

#print D
#print rows
#print cols


sys.stderr.write("Determining subset matrix\n")
E = {}

for r in rows:
  #print r
  E[r] = {}
  for c in cols:
    #print "  ", c
    E[r][c] = D[r][c]
    #print D[r][c]

#print_matrix(E)  # This does not preserve row/col order
sys.stderr.write("Writing out result\n")
for c in cols:
  sys.stdout.write("\t%s" % c)
sys.stdout.write("\n")

for r in rows:
  if OMIT_ZERO_ROWS:
    sum = 0
    for c in cols:
      sum += int(E[r][c])
    if sum == 0:
      continue
  # If sum > 0, we'll still be in this iteration of the loop. Below will be executed.
  sys.stdout.write(r)
  for c in cols:
    sys.stdout.write("\t%s" % E[r][c])
  sys.stdout.write("\n")
