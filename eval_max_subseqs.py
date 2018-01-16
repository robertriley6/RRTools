# Find maximum scoring intervals of a multiple alignment where each position is assigned some score.
 
# Read in a multiple alignment
# For each column, assign a score
# From this you get a sequence of scores.
# Use the method of Tompa and Ruzzo (1999) to find all maximum scoring subsequences

import sys
import re
import os
from seq_functions import *

try:
  print "reading alignment..."
  aln = read_fasta_alignment(sys.argv[1])
  print "read alignment %s" % sys.argv[1]
except:
  print "Supply a valid FASTA alignment file"
  sys.exit()

# Minimum length of conserved region we will consider
MIN_LENGTH = 1


# Define functions
def process_cmd(path):
  a = re.split("/", path)
  filename = a[-1]
  b = re.split("\.", filename)
  prefix = b[0]
  return prefix
  

# PROCESS INPUT
prefix = process_cmd(sys.argv[1])  # Return a prefix to use as output file name

# The length of alignment
N = aln.length()

# The sequence of positionwise conservation scores
#print "Computing position scores..."
scores = []
cumulative_scores = []

print "Calculating positional scores..."
for i in range(N):
  col = aln.get_column(i + 1)
  s = position_conservation(col)
  scores.append(s)

# Test
#scores = [4, -5, 3, -3, 1, 2, -2, 2, -2, 1, 5]
print "normalizing..."
scores = [scores[i] - .5 for i in range(len(scores))]

#for s in scores:
#  print s
#sys.exit()

print "calculating cumulative scores"
cumulative_scores = [sum(scores[0:i+1]) for i in range(len(scores))]

I = eval_max_subseq(scores)

# All subsequences in the list are now maximal.  Output them.
print "writing results"
result_file = open(prefix + '.cons', 'w')
print_results(I, MIN_LENGTH, prefix, result_file)

print "writing scores"
score_file = open(prefix + '.score', 'w')
print_scores(scores, cumulative_scores, score_file)

if not result_file.tell():
  os.remove(result_file.name)
  os.remove(score_file.name)
  print "No results for %s" % result_file.name

result_file.close()
score_file.close()
