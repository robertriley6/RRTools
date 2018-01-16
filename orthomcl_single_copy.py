# Extract single-copy clusters from OrthoMCL run, e.g. file 'groups.txt'
#
# File has one cluster per line, lines look like:

# cluster16469: asp|1176580 cryp|42756 sta|8442 tri|75434
# cluster16470: sch|2664171 dic|139591 fomp|1166292 glo|128976
# cluster16471: asp|1177110 cryp|357688 sta|3771 tri|22774

# Assume your proteins are named e.g. 'asp|1176580' where 'asp' is a species 
# name and '1176580' is a unique identifier for that protein in asp, etc.

# Write out those clusters that have one and only one protein from each species

import sys
import re

TRUE = 1
FALSE = 0

try:
  f = open(sys.argv[1], 'r')
except:
  print "Supply a file, e.g. 'groups.txt' from an OrthoMCL run"
  sys.exit()

clusters = {}

# Read in the clusters
for line in f.readlines():
  m = re.match("^(\w+):", line)
  if m:
    cluster_id = m.group(1)
    clusters[cluster_id] = re.findall('(\w+)\|(\d+)', line)

# Determine what species are in the clusters
species = set([s for c in clusters.itervalues() for s, p in c])


# Find the cluster, where one and only one protein is contributed by each species
single_copy_clusters = {}
N = len(species)
for k, c in clusters.iteritems():
  if len(c) == N:
    species_in_cluster = set([s for s, p in c])  # Give the list of species in a cluster
    if len(species_in_cluster) == N:
      single_copy_clusters[k] = c

sys.stderr.write("%d single-copy clusters obtained\n" % len(single_copy_clusters))

for k, c in single_copy_clusters.iteritems():
  sys.stdout.write("%s:" % k)  # Write the cluster id
  for i in range(N):  # Write each protein in the cluster
    s, p = c[i]
    sys.stdout.write(" %s|%s" % (s, p))
  sys.stdout.write("\n")
