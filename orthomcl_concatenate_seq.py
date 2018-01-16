# Read in a list of single-copy clusters from some MCL-based clustering method 
# (e.g. OrthoMCL), and concatenate the sequences for each organism.  For each
# organism, write out one FASTA sequence containing the concatenated sequences,
# in the same cluster order for each organism.

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

all_fastas = {}
clusters = []

try:
  try:
    fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
    for fasta in fasta_sequences:
      name, sequence = fasta.id, str(fasta.seq)
      all_fastas[name] = fasta  # fasta is a SeqRecord object.  See http://biopython.org
  except:
    sys.stderr.write("ERROR: supply a valid FASTA file\n")
    sys.exit()
  try:
    f = open(sys.argv[2], 'r')
    for line in f.readlines():
      c = re.findall('(\w+)\|(\d+)', line)
      if c:
        clusters.append(c)
  except:
    sys.stderr.write("ERROR: supply a cluster file\n")
    sys.exit()
  try:
    out_file = open(sys.argv[3], 'w')
  except:
    sys.stderr.write("ERROR opening %s\n" % sys.argv[3])
    sys.exit() 
except:
  print "ARGUMENTS: <input fasta file> <single copy clusters file> <output fasta file name>"
  sys.exit()
 
all_species = set([s for c in clusters for s, p in c])

concat_fastas = []  # Initialize dictionary

for species in all_species:
  concat_seq = ''
  for c in clusters:
    for s, p in c:
      if s == species:
        concat_seq += str(all_fastas[s + '|' + p].seq)
  concat_fastas.append(SeqRecord(Seq(concat_seq, IUPAC.protein), id=species, description=""))

SeqIO.write(concat_fastas, out_file, "fasta")
