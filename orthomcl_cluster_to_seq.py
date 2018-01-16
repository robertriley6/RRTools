# Read in a list of single-copy clusters from some MCL-based clustering method 
# (e.g. OrthoMCL), and concatenate the sequences for each organism.  For each
# such cluster, create an output file.  In each of those files (clusters), for
# each organism, write out one FASTA sequence containing the concatenated 
# sequences, in the same cluster order for each organism.

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

all_fastas = {}
clusters = {}

try:
  try:
    sys.stderr.write("Reading fastas...\n")
    fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
    for fasta in fasta_sequences:
      name, sequence = fasta.id, str(fasta.seq)
      all_fastas[name] = fasta  # fasta is a SeqRecord object.  See http://biopython.org
  except:
    sys.stderr.write("ERROR: supply a valid FASTA file\n")
    sys.exit()
  try:
    sys.stderr.write("Reading clusters...\n")
    f = open(sys.argv[2], 'r')
    for line in f.readlines():
      m = re.search("^(\w+):", line)
      if m:
        k = m.group(1)  # The cluster id
        c = re.findall('(\w+)\|(\d+)', line)
        if c:
          clusters[k] = c
  except:
    sys.stderr.write("ERROR: supply a cluster file\n")
    sys.exit()
except:
  print "ARGUMENTS: <input fasta file> <single copy clusters file>"
  sys.exit()
 
all_species = set([s for c in clusters.itervalues() for s, p in c])


for k in clusters.keys():  # For each cluster
  out_file = open(k + ".fasta", 'w')
  out_fastas = []
  for s, p in clusters[k]:
    name = s + "|" + p  # Should reconstruct name of fasta header
    out_fastas.append(all_fastas[name])
  SeqIO.write(out_fastas, out_file, "fasta")
  out_file.close()
