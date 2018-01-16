import sys
from re import compile, I, search, split, sub
from stats import *
from math import *
from MySQLdb import *
from httplib import *
from classes import *

FALSE = 0
TRUE = 1


def read_list(path):
  # Read a generic list from a file
  f = open(path, 'r')
  L = []
  for line in f.readlines():
    pat = search(r'^(.+)$', line)  # Notice that list elements can have whitespace
    if pat:
      L.append(pat.group(1))
  if len(L) == 0:
    print "WARNING: list empty"
  return L



def write_list(L):
  for l in L:
    sys.stdout.write(l + "\n")



def read_tab_delimited(path):
  f = open(path, 'r')
  data = []
  for line in f.readlines():
    data.append(split("\t", line[:-1]))
  return data



def read_fasta(path):  # returns dict
  seq = 0  # for TRUE/FALSE logic purposes
  f = open(path)
  data = {}
  for line in f.readlines():
    h_pat = search(r'^>(.+)$', line)
    d_pat = search(r'^([^>]+)$', line)

    if h_pat:
      if seq:  # only not true when first header is read
        data[seq[0]] = seq[1]
      seq = [h_pat.group(1), '']

    elif d_pat:
      seq[1] += d_pat.group(1).replace('\n', '')

  # Process the last sequence
  data[seq[0]] = seq[1]

  f.close()
  return data



def write_fasta(S, n):  # S == dict of sequences, n == number of character per line
  for s in S.keys():
    header = '>' + s
    seq = S[s]
    sys.stdout.write(header + "\n")
    L = len(seq)
    i = 0  # you have to initialize i to some value in case below ratio == 0
    while i <= L/n:
      sys.stdout.write(seq[i*n:(i+1)*n] + "\n")
      i += 1



def read_matrix(path):  # Make a dictionary
  # D[row][col]
  # Input looks like
  #                 pc      pp      sc      sl      an
  # Glyco_hydro_61  38      33      66      14      46
  # Glycos_transf_2 26      44      25      23      26

  # Assumes tab-delimited

  # Should also handle empty values "\t\t"

  # We want to end up with D = {'Glyco_hydro_61':{'pc':38, 'pp':33, ...}, ...}
  f = open(path, 'r')

  D = {}  # row_id -> [0, 1, 55,.., N]
  COLS_KNOWN = FALSE  # We don't yet know the names of the columns
  DATA_SEEN = TRUE  # Whether or not we've read any data lines
  for line in f.readlines():
    col_pat = search(r'^\t(.+)$', line)  # Get the first line's column ids, whatever they might be.  Tab at beginning.
    dat_pat = search(r'^(\S[\S ]+)\t(.+)$', line)

    if col_pat:
      cols = split('\t', col_pat.group(1))
      L = len(cols)
      COLS_KNOWN = TRUE

    if dat_pat and COLS_KNOWN:
      row = dat_pat.group(1)
      cts = split('\t', dat_pat.group(2))
      if len(cts) == L:
        D[row] = {}
        for c in cols:
          i = cols.index(c)
          D[row][c] = cts[i]
      else:
        print "ERROR: vector of strange length"
        sys.exit()
  if not COLS_KNOWN:
    print "ERROR: header row suspect"
  if not DATA_SEEN:
    print "ERROR: data rows suspect"
  return D



def clean_list(L):
  for l in L:
    if l == '':
      L.remove(l)
  return L



def list_to_int(L):
  for i in range(len(L)):
    L[i] = int(L[i])
  return L


# Read in pairs
# a0 b0
# a1 b1
# to a list
def read_pairs(path):
  L = read_list(path)
  M = [split("\s+", l) for l in L]
  return M



# make dict a0 -> b0 etc.
def read_pairs_to_dict(path):
  L = read_list(path)
  D = {}
  for l in L:
    a, b = split("\s+", l)
    D[a] = b
  return D



# same but b0 etc. becomes int
def read_pairs_to_dict_int(path):
  L = read_list(path)
  D = {}
  for l in L:
    a, b = split("\s+", l)
    D[a] = int(b)
  return D



# Make dict of (a, b) -> v mappings
def read_valued_pairs_to_dict(path):
  L = read_list(path)
  D = {}
  for l in L:
    a, b, v = split("\s+", l)
    # Intitialize dicts as appropriate
    if not D.has_key(a):
        D[a] = {}
    # Print warning if we overwrite an existing D[a][b]
    if D[a].has_key(b):
      print "WARNING: duplicate pair: %s %s" % (a, b)
    # Symmetrically assign v
    D[a][b] = v
  return D
  


# Turn vector to z-scores
def Zify_vector(X):
  L = len(X)
  
  Z = []  # The vector to put indicators in
  for i in range(L):
    A = [X[j] for j in range(L) if j != i]  # Leaves X[i] out of distribution for determining Z
    Z.append(z_score(X[i], A))
  return Z



def binarize_vector(X):
  A = []
  for x in X:
    if x == 0:
      A.append(0)
    else:
      A.append(1)
  return A



def binarize_by_z(X, z):
  B = []
  for x in X:
    if x > z:
      B.append(1)
    else:
      B.append(0)
  return B



def pfam2go_dict(path):
  # The lines we seek look like this:
  # Pfam:PF00001 7tm_1 > GO:G-protein coupled receptor protein signaling pathway ; GO:0007186
  # Pfam:PF00006 ATP-synt_ab > GO:hydrolase activity, acting on acid anhydrides, catalyzing transmembrane movement of substances ; GO:0016820
  # etc.
  
  f = open(path, 'r')

  D = {}

  for line in f.readlines():
    pat = search(r'^Pfam:PF\d+\s(\S+) > GO:(.+)\s;\s(GO:\d+)$', line)
    if pat:
      go_id = pat.group(3)
      go_desc = pat.group(2)
      pf_id = pat.group(1)
      if not D.has_key(pf_id):
        D[pf_id] = [(go_id, go_desc)]
      else:
        D[pf_id].append((go_id, go_desc))
    else:
      pass
  return D



# Invert dict of dicts
# Usually in the case of a square matrix represented as a dictionary, we want to transpose
def invert_dict(D):
  E = {}
  X = []
  Y = []

  for x in D.keys():
    if x not in X:
      X.append(x)
    for y in D[x].keys():
      if y not in Y:
        Y.append(y)

  for y in Y:
    if not E.has_key(y):
      E[y] = {}
    for x in X:
      try:
        E[y][x] = D[x][y] 
      except:
        pass
  return E  



# Print a matrix to stdout
# Format:
#
# 	Ab	Cc	Cn
# GH1	0	2	5
# GH2	9	7	3

# Below requires all nested dicts have the same set of keys.
def print_matrix(M):  # Takes a dict of dicts as argument
  rows = M.keys()
  cols = M[rows[0]].keys()
  rows.sort() # Better sorted than random order
  cols.sort()
  
  # First line
  for c in cols:
    sys.stdout.write('\t%s' % c)
  sys.stdout.write('\n')

  # Data lines
  for r in rows:
    sys.stdout.write('%s' % r)
    for c in cols:
      sys.stdout.write('\t%s' % M[r][c])
    sys.stdout.write('\n')
    sys.stdout.flush()


# Same as print_matrix, but no requirement that all nested dicts have same set of keys
def print_sparse_matrix(M):  # Takes a dict of dicts as argument
  rows = M.keys()
  cols = []
  for r in rows:
    for c in M[r].keys():
      if c not in cols:
        cols.append(c)
  rows.sort() # Better sorted than random order
  cols.sort()
  
  # First line
  for c in cols:
    sys.stdout.write('\t%s' % c)
  sys.stdout.write('\n')

  # Data lines
  for r in rows:
    sys.stdout.write('%s' % r)
    for c in cols:
      try:
        sys.stdout.write('\t%s' % M[r][c])
      except:
        sys.stdout.write('\t0')
    sys.stdout.write('\n')
    sys.stdout.flush()



def reverse_complement(seq):
  reverse = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
  complement = ''
  for res in seq:
    complement += reverse[res]
  reverse_complement = complement[::-1]
  return reverse_complement



def log_transform_dict(X):
  Y = {}
  for a in X.keys():
    Y[a] = {}
    for b in X[a].keys():
      Y[a][b] = log10(float(X[a][b]))  # in case we're handed a data that isn't floats
  return Y



# KEGG functions
def get_kegg_object(obj_type, name):
  # curl -s http://rest.kegg.jp/get/cpd:C00042
  # compound = cpd e.g. C00042
  # pathway = path e.g. map00020
  # enzyme = ec e.g. 1.2.1.16
  OBJ_FOUND = FALSE
  conn = HTTPConnection("rest.kegg.jp")
  conn.request("GET", "/get/" + obj_type + ":" + name)
  r1 = conn.getresponse()
  data = r1.read()
  lines = split("\n", data)
  for line in lines:
    pat = search(r'^NAME\s+(.+)$', line)
    if pat:
      value = pat.group(1)
      OBJ_FOUND = TRUE
      break

  if not OBJ_FOUND:
    sys.exit()

  return value
 


def get_kegg_objects(obj_type, names):
  # curl -s http://rest.kegg.jp/get/cpd:C00042
  # compound = cpd e.g. C00042
  # pathway = path e.g. map00020
  # enzyme = ec e.g. 1.2.1.16
  result = []

  conn = HTTPConnection("rest.kegg.jp")
  for name in names:
    OBJ_FOUND = FALSE
    conn.request("GET", "/get/" + obj_type + ":" + name)
    r1 = conn.getresponse()
    data = r1.read()
    lines = split("\n", data)
    for line in lines:
      pat = search(r'^NAME\s+(.+)$', line)
      if pat:
        value = pat.group(1)
        OBJ_FOUND = TRUE
        result.append(value)
        break

    if not OBJ_FOUND:
      sys.exit()

  return result
 


# read enzyme.dat filed downloaded from ENZYME and get dict with ec -> description mappings, e.g.
# 1.1.1.9 -> D-xylulose reductase.
# typical data lines 
# ID   1.1.1.1
# DE   Alcohol dehydrogenase.
# //

def get_ec_desc(path):
  D = {}  # ec -> desc
  ID = ''  # value from ID field in enzyme.dat
  DE = ''  # value from DE field in enzyme.dat

  f = open(path, 'r')

  for line in f.readlines():
    id_pat = search(r'^ID\s+(\d+\.\d+\.\d+\.\d+)$', line)
    de_pat = search(r'^DE\s+(\S.+)$', line)
    end_pat = search(r'^\/\/$', line)
    
    if id_pat:
      ID = id_pat.group(1)

    if de_pat:
      DE = de_pat.group(1)

    if end_pat:
      # process data if appropriate
      if ID != '' and DE != '':  # if it's after the first time we see '//', because there is one before the first record
        D[ID] = DE
  return D



###################
# GENERAL FUNCTIONS
# FOR THIS AND THAT
###################

# Do the intervals overlap at all?
def overlaps(M, N):
  a, b = M
  c, d = N
  return (min(a,b) <= c or min(a,b) <= d) and (max(a,b) >= c or max(a,b) >= d)

def contains(M, N):
  a, b = M
  c, d = N
  return (min(a,b) <= min(c,d)) and (max(a,b) >= max(c,d)) 

def is_contained_in(M, N):
  a, b = M
  c, d = N
  return (min(a,b) >= min(c,d)) and (max(a,b) <= max(c,d)) 

