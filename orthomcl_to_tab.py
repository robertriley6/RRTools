# Take groups.txt file from OrthoMCL run and make input file for Count, etc.
# Cluster file lines look like
# cluster15017: ser|1209981 tra|71397 tra|46569 dac|111206

# Count input data must look like:

# family Aerpe Arcfu Calma Censy Halma Halsp
# arCOG00001 1 2 1 

import sys
import re

try:
  f = open(sys.argv[1], 'r')
except:
  print "Supply a groups.txt file from an OrthoMCL run"
  sys.exit()

clusters = {}
C = []  # List of cluster ids to preserve order we read in clusters

# Read in the clusters
for line in f.readlines():
  m = re.match("^(\w+):", line)
  if m:
    cluster_id = m.group(1)
    clusters[cluster_id] = re.findall('(\w+)\|(\d+)', line)
    C.append(cluster_id)
# Determine what species are in the clusters

species = set([s for c in clusters.itervalues() for s, p in c])

# Write out the header line
sys.stdout.write("family")
for s in species:
  sys.stdout.write("\t%s" % s)
sys.stdout.write("\n")

# Write out the data
for i in C:
  sys.stdout.write("%s" % i)
  for s in species:
    n = len([org for org, prot in clusters[i] if org == s])
    sys.stdout.write("\t%d" % n)
  sys.stdout.write("\n")
