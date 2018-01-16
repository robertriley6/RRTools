# Read in a list of gene family annotations in several genomes, e.g.
# 113829|Cersu1	GH63
# 348313|Cerun2	GH6
# 391571|Cerun2	GH63
# 1159900|Clibor1	GH6
# 1180007|Clibor1	GH63
# 747039|Daequ1	GH63
# 152125|Dicsq1	GH6

# Read in a list of organisms making up a meaningful (e.g. brown rot)

# subest of the organisms in the matrix, e.g.

  # Cersu1
  # Dicsq1

# Tabulate for each family
#
#   a - has pfam, in subset
#   b - has pfam, in background distribution
#   c - doesn't have pfam, in subset
#   d - doesn't have pfam, in background distribution

# Contingency table for each family is
#         Atlantic  Indian
#        whales     8        2
#        sharks     1        5

# a, b, c, and d are all mutually exclusive.


import sys
import scipy.stats as stats
from functions import read_pairs, read_list
from re import split


if len(sys.argv) < 4:
  print "ARGS: <annotation file> <subset org list file> <type of test: 0 = depletion, 1 = enrichment, 2 = two-tailed> <optional significance threshold>"
  sys.exit()

annotation_path = sys.argv[1]
subset_path = sys.argv[2]
alt_arg = sys.argv[3]
try:
  p_threshold = float(sys.argv[4])
except:
  p_threshold = 0.05

# Read list of subset organisms
P = set(read_list(subset_path))

#Q = set(read_list(background_path)).difference(P)

# What kind of test do you want to perform
if alt_arg == '0':
  alt = 'less'
elif alt_arg == '1':
  alt = 'greater'
elif alt_arg == '2':
  alt = 'two-sided'
else:
  alt = 'two-sided'  # A reasonable default


class Family:
  def __init__(self, name):
    self.name = name
    self.a = 0  # Number of proteins with this annotation, in subset
    self.b = 0  # Number of proteins with this annotation, not in subset
    self.c = 0  # Number of proteins without this annotation, in subset
    self.d = 0  # Number of proteins without this annotation, not in subset
    self.protlist = []  # List of proteins with this annotation


F = {}  # Families. family_id -> Family object
O = {}  # org -> protein -> family annotations

family_list = []

# Read annotations
sys.stderr.write("Reading in annotations from %s\n" % annotation_path)
for prot, family in read_pairs(annotation_path):  # 12345|Copci -> GH6
  pid, org = split("\|", prot) 
  #print org, pid, family
  if org not in O.keys():
    O[org] = {}
  if pid not in O[org]:
    O[org][pid] = []  # List of domains
  O[org][pid].append(family)
  if not F.has_key(family):
    F[family] = Family(family)  # Initialize

# For each family
sys.stderr.write("Tabulating annotations...\n")
for f in F.keys():
  for org in O.keys():
    for prot in O[org].keys():
      if f in O[org][prot] and org in P:
        F[f].a += 1
      elif f in O[org][prot] and org not in P:
        F[f].b += 1
      elif f not in O[org][prot] and org in P:
        F[f].c += 1
      elif f not in O[org][prot] and org not in P:
        F[f].d += 1
      else:   # Should never happen
        print "ERROR:", f, org, prot
        sys.exit()


alpha = len(F.keys())  # Bonferroni correction

sys.stderr.write("Doing Fisher's test\n")
sys.stdout.write("FAMILY\tA\tB\tC\tD\tRATIO\tPVAL\n")
for f in F.keys():
  a = F[f].a
  b = F[f].b
  c = F[f].c
  d = F[f].d

  oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]], alt)
  corrected_p = min(alpha * pvalue, 1.0)
  if corrected_p <= p_threshold:
    sys.stdout.write("%s\t%d\t%d\t%d\t%d\t%f\t%f\n" % (F[f].name, a, b, c, d, oddsratio, corrected_p))
    sys.stdout.flush()

