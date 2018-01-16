####################################################
# VARIOUS FUNCTIONS FOR DEALING WITH BIOLOGICAL DATA
####################################################

import sys
import re
from re import findall, search, split
from math import log
from seq_classes import *
from pickle import load


###############
# THIS AND THAT
###############

TRUE = 1
FALSE = 0
ERROR = -1


###########
# FUNCTIONS
###########


# What percentage of an alignment is actual residues (not gaps)?
def alignment_density(aln):
  n_res = 0
  n_char = 0
  for s in aln.sequences:
    n_res += len(findall('[a-zA-Z]', s.data))
    n_char += len(s.data)
  return float(n_res) / n_char  


# Read a FASTA format multiple alignment.  Return an Alignment object
# Check first to see that it appears to be an alignment (all sequences
# same length including gaps)
def read_fasta_alignment(file_name):
  f = open(file_name, 'r')

  aln = Alignment()

  header = ''
  id = ''
  seq = ''
  
  for line in f.readlines():
    header_pat = search(r'^>(.+\S)\s*$', line)
    data_pat = search(r'^\s*([a-zA-Z\-]+)\s*$', line)
    
    if header_pat:
      # Suppose it's not the first one we've seen
      # Dump the current header, seq into a Seq object
      if len(seq) > 0:  # Meaning you probably have just read the first line
        S = Seq(header, seq)   
        S.id = id 
        aln.sequences.append(S)
      # and reset	
      header = header_pat.group(1)
      id_pat = search('^(\w+).*', header)
      if id_pat:
        id = id_pat.group(1)
      print header, id
      seq = ''	

    if data_pat:
      seq += data_pat.group(1).upper()
  # Loop ends, store the last one
  S = Seq(header, seq)   
  S.id = id 
  aln.sequences.append(S)

  # And now a reality check
  L = len(aln.sequences[0].data)
  for s in aln.sequences:
    k = len(s.data)
    if L != k:
      print "ERROR: different length sequences in alignment."
      print "This doesn't look like an alignment."
      sys.exit()
  if len(aln.sequences) == 0:
    print "WARNING: no sequences in alignment"
  # Make sure all is in upper case
  return aln



# Read one or more fasta formatted sequences from a file
def read_fasta(file_name):
  f = open(file_name, 'r')
  
  seq_list = []
  
  header = ''
  seq = ''
  
  for line in f.readlines():
    header_pat = search(r'^>(.+)$', line)
    data_pat = search(r'^([a-zA-Z\-]+)$', line)
    
    if header_pat:
      # Suppose it's not the first one we've seen
      # Dump the current header, seq into a Seq object
      if len(seq) > 0:  # Meaning you probably have just read the first line
        seq_list.append(Seq(header, seq))
      # and reset	
      header = header_pat.group(1)
      seq = ''	

    if data_pat:
      seq += data_pat.group(1)
  # Loop ends, store the last one
  seq_list.append(Seq(header, seq))    

  return seq_list


# What is the alignment length?
def alignment_length(aln):
  L = len(aln.sequences[0].data)
  for s in aln.sequences:
    k = len(s.data)
    if L != k:
      print "ERROR: different length sequences in alignment"
      sys.exit()
  return L    


# Get a subalignment, Xab, from alignment X spanning the interval (a,b)  
def get_subalignment(X, a, b):

  # Check for errors
  if b < a:
    print "ERROR: b shouldn't be less than a"
    return ERROR

  elif a < 0:
    print "ERROR: indices must be positive"
    return ERROR
    
  else:
    # Continue
    Xab = Alignment()
    for s in X.sequences:

      new_header = update_interval(s, a, b)
      new_seq = s.data[a-1:b]   # This is how you get a,b inclusive, where 1 (not 0) is the number we start from
      Xab.sequences.append(Seq(new_header, new_seq))
    return Xab  



#  With a sequence S  from an alignment looking maybe like this '-----AACCC---TT---'
#  and knowing that it corresponds to coordinates x,y on a genome
#  judging from the header looking like
#  >MA104 7000000046786761:   11170-   13699 +
#  come up with the new set of genomic coordinates represented in the subalignment
#  X(a,b)
def update_interval(S, a, b):
  h_pat = search(r'(^\S+)\s+.*:\s+(\d+)-\s*(\d+)\s+([+\-])', S.header)
  if h_pat:
    genome_id = h_pat.group(1)
    x = int(h_pat.group(2))
    y = int(h_pat.group(3))
    strand = h_pat.group(4)
    m = x + (a - 1) - num_gaps(S.data[0:a-1])
    n = m + (b - a) - num_gaps(S.data[a-1:b])
    new_header = "%s : %d- %d %s"  % (genome_id, m, n, strand)
  else:
    print "unrecognizable header: %s" % S.header
    sys.exit()
  return new_header    



def translate_interval(S, x, y, a, b):
  m = x + (a - 1) - num_gaps(S.data[0:a-1])
  n = m + (b - a) - num_gaps(S.data[a-1:b])
  return (m,n)

# How many gaps in a string
def num_gaps(seq_str):
  n = 0
  for s in seq_str:
    if s == '-':
      n += 1
  return n



# Read a ptt format file as in ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Mycobacterium_tuberculosis_H37Rv/NC_000962.ptt
def read_ptt_annot(file):
  annotations = []  # Ordered
  f = open(file, 'r')
  for line in f.readlines():
    pat = search(r'^(\d+)\.\.(\d+)\t([\+\-]{,1})\t(\d+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)$', line)
    if pat:
      s = int(pat.group(1)) 
      e = int(pat.group(2)) 
      if pat.group(3):  # If strand info was given in ptt file, usual, but not always
        strand = pat.group(3)
      else:
        strand = ''
      length = int(pat.group(4))
      pid = pat.group(5)
      gene = pat.group(6)
      syn = pat.group(7)
      code = pat.group(8)
      cog = pat.group(9)
      prod = pat.group(10)
      annotations.append(PTT_annot(s, e, strand, length, pid, gene, syn, code, cog, prod))
    else:
      pass
  return annotations



def read_gff(filename):
  G = []
  f = open(filename, 'r')
  for line in f.readlines():
    pat = search('^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(.*)\t(.*)\t(.*)\t(.+)$', line)
    if pat:
      seqname = pat.group(1)
      source = pat.group(2)
      feature = pat.group(3)
      start = int(pat.group(4))
      end = int(pat.group(5))
      score = pat.group(6)
      strand = pat.group(7)
      frame = pat.group(8)
      attributes = pat.group(9)
      if '#' in attributes:
        comments = search(r'.*(#.*)$', attributes).group(1)
      else:
        comments = ''
      G.append(GFF_annot(seqname, source, feature, start, end, score, strand, frame, attributes, comments))
  return G
 


def make_annot_dict(annotations):
  annot_dict = {}
  for n in annotations:
    a = n.start
    b = n.end
    annot_dict[(a,b)] = n
  return annot_dict



# How many nucleotides are in a sequence, regardless of gaps 
def n_bases(seq):
  n_nuc = 0
  for i in seq.data:
    if i != '-' and i != ' ':  # If it's not a gap or whitespace
      n_nuc += 1
  return n_nuc



def n_bases_gaps(seq):
  n_nuc = 0
  n_gap = 0
  for i in seq.data:
    if i != '-' and i != ' ':  # If it's not a gap or whitespace
      n_nuc += 1
    else:
      n_gap += 1
  return n_nuc, n_gap



# Does an aligment column (provided as a string) have any gaps?
def has_gap(col):
  for i in col:
    if i == '-':
      return TRUE
  return FALSE



# Is an alignment column conserved in all sequences?
def is_conserved(col):
  a = col[0]  # The first character, arbitrary.
  for i in col[1:]:
    if i.lower() != a.lower():  # Correct for case
      return FALSE
  return TRUE



def string_entropy(X):
  base = 2
  counts = {}
  for x in X:
    if not counts.has_key(x):
      counts[x] = 1
    else:
      counts[x] += 1
  n = len(X)
  sum = 0.0
  for x in X:
    px = counts[x] / float(n)         
    sum += px * log(px, base)
  return -1 * sum



def position_conservation(col):
  N = len(col)
  denom = N * (N - 1) / 2
  numer = 0
  for i in range(N):
    for j in range(N):
      if i > j:
        if col[i] == col[j] and not (col[i] == '-' or col[j] == '-'):
          numer += 1
  return numer / float(denom)



# Look at an alignment column and infer conservation based on position identity,
# weighted by the distance from an alignment-derived matrix
def position_phylo_cons(col, D, genomes, ref_gen):
  N = len(col)  # Keeping in mind that the ref sequence is not an informant
  # The denominator
  w_sum = 0.0
  for g in genomes:
    w_sum += D[id_mappings[ref_gen]][id_mappings[g]]
  ref_base = col[genomes.index(ref_gen)]
  if ref_base == '-':  # If it's a gap
    return 0.0
  # Otherwise, calculate the numerator 
  sum = 0.0
  for i in range(N):
    if genomes[i] != ref_gen:
      if col[i] == ref_base:
        sum += D[id_mappings[ref_gen]][id_mappings[genomes[i]]]
  return sum / w_sum         



def dictify_file_list(f_list):
  D = {}
  for f in f_list:
    arr = split("\.", f)
    if len(arr) == 2:
      a = arr[0]
      b = arr[1]
      D[a] = b
    else:
      print "WARNING: strange looking file: %s" % f
  return D 



def parse_distance_matrix(file_name):
  # Parse a ClustalW Phylip-format distance matrix
  # Put numbers in a Numeric data structure
  DIMENSION_KNOWN = FALSE  # The file has an integer on first line telling us dimension of matrix
  FIRST_ROW_SEEN = FALSE
  f = open(file_name, 'r')
  i = 0
  for line in f.readlines():
    if not DIMENSION_KNOWN:  # Then we should look for it
      header_pat = search(r'^\s+(\d+)\s+$', line)
      if header_pat:
        dim = int(header_pat.group(1))
        D = zeros((dim, dim), float)
        DIMENSION_KNOWN = TRUE
    else:
      data_start_pat = search(r'^(\S+\s+.+\s)$', line)
      data_more_pat = search(r'^(\s+.+\s)$', line)
      if data_start_pat:
        if not FIRST_ROW_SEEN:  # If this is the first row 
          row = split('\s+', data_start_pat.group(1))[1:-1]
          FIRST_ROW_SEEN = TRUE
        else:
          # This means it is the start of a new row, but not the first we've read
          # Put old row's data in matrix
          for j in range(len(row)):
            D[i][j] = float(row[j])
          i +=1  # Only increment row counter after we've loaded previous row
        # Start the new row     
        row = split('\s+', data_start_pat.group(1))[1:-1]
      elif data_more_pat:
        row += split('\s+', data_more_pat.group(1))[1:-1]
  # Write the last one
  for j in range(len(row)):
    D[i][j] = float(row[j])
  return D



def get_distance_dict(file_name):
  from re import split
  # Parse a ClustalW Phylip-format distance matrix
  # Put numbers in a dictionary keyed by pairs of protein ids 
  DIMENSION_KNOWN = FALSE  # The file has an integer on first line telling us dimension of matrix
  f = open(file_name, 'r')
  row = 0  # A sentinel to indicate whether we have the first row or not later
  data = []  # Store the data on the lines
  key_list = []  # Ordered keys. Think of as the values on i and j axes
  for line in f.readlines():
    if not DIMENSION_KNOWN:  # Then we should look for it
      header_pat = search(r'^\s+(\d+)$', line)
      if header_pat:
        dim = int(header_pat.group(1))
        DIMENSION_KNOWN = TRUE
        print "%d dimensional distance matrix" % dim
    else:
      data_start_pat = search(r'^(\S+)\s+(.+\s)$', line)
      data_more_pat = search(r'^(\s+.+\s)$', line)
      if data_start_pat:
        # If there's an old row, store it
        if row:
          key_list.append(id)
          data.append(row)
        # Start the new row     
        id = data_start_pat.group(1)
        row = split('\s+', data_start_pat.group(2))[0:-1]
      elif data_more_pat:
        row += split('\s+', data_more_pat.group(1))[1:-1]
  # store the last one
  key_list.append(id)
  data.append(row)
  #Put it all in a dict
  D = {}
  for i in key_list:
    x = key_list.index(i)
    D[i] = {}
    for j in key_list:
      y = key_list.index(j)
      D[i][j] = float(data[x][y])
  return D



# Functions for implementing Tompa and Ruzzo's algorithm for detecting maximum scoring subsequences

def compute_L(scores, k):
  # Initialize L 
  L = 0
  
  # Compute L 
  # ... cumulative total L_j of all scores up to but not including the leftmost score of I_j
  i = 0
  while i < k:  # Recall that I = (left, right, score)
    L += scores[i]
    i += 1
  return L



def compute_R(scores, k):
  # initialize R
  R = 0

  # Compute R
  # ... and the total R_j up to and including the rightmost score of I_j
  i = 0
  while i <= k:
    R += scores[i]
    i += 1
  return R



def compute_subseq_score(scores, l, r):
  s = 0
  i = l
  while i <= r:
    s += scores[i]
    i += 1
  return s



def print_scores(S, C, file):
  for i in range(len(S)):
    file.write("%d\t%.2f\t%.2f\n" % (i, S[i], C[i]))



def print_results(I, MIN_LENGTH, label, file):
  i = 0
  for I in I:
    if (I.right +1) - I.left >= MIN_LENGTH:
      file.write("%d\t%d\t%s %.2f\n" % (I.left + 1, I.right + 1, label + '_' + str(i), I.score))
      i += 1



def eval_max_subseq(scores):
  N = len(scores)
  # RUN THE ALGORITHM
  print "Running algorithm..."
  # Maintain cumulative total of scores read so far
  # Maintain ordered list of disjoint subsequences
  I = []  # The list is initially empty
  
  # Read scores from left
  k = 0  
  while k < N:
    if scores[k] <= 0.0:
      pass  # A nonpositive score requires no special processing when read
      k += 1
    else:
      # Incorporate positive score into new subsequence I_k with length one, 
      try: I_k  # Test for the existence of I_k.  Are we inheriting it from the last loop? Or do we initialize a new one of length 1?
      except: I_k = Subseq(k, k, scores[k])  # (left, right, score)
      # Either way now, we have an I_k
      # and integrate into list as follows
  
      # 1. The list is searched from right to left for the maximum value of j
      #    satisfying L_j < L_k
      L_k = compute_L(scores, I_k.left)
      for j in range(len(I)-1, -2, -1):  # Count down from len(I), basically.  Notice our last j will be -1, to signify no suitable j found, if loop not broken
        if j < 0:
          break
        L_j = compute_L(scores, I[j].left)
        if L_j < L_k:
          break  # Stop, you found the leftmost j that satisfied the condition.  j will remain set to the appropriate value  
  
      # If there is no such j, then add I_k to the end of the list
      if j < 0:
        I.append(I_k) 
        del I_k
        k += 1
      else:
        # If there is such a j, and R_j >= R_k, then add I_k to the end of the list
        R_j = compute_R(scores, I[j].right)
        R_k = compute_R(scores, I_k.right)
        if R_j >= R_k:
          I.append(I_k)
          del I_k
          k += 1
  
        # Otherwise (there is such a j, but R_j < R_k), extend subsequence I_k to the left
        # to encompass everything up to and including the leftmost score in I_j.  Delete
        # subsequences I_j, I_{j+1}, ..., I_{k-1} from the list (none of them is maximal)
        # and reconsider the newly extended subsequence I_k (now renumbered I_j) as in step 1.
        else:
          I_k.left = I[j].left
          I_k.score = compute_subseq_score(scores, I_k.left, I_k.right) # Update the score to reflect new bounds
          del I[j:k]  # This does not include k, stops at k-1
  # All subsequences in the list are now maximal.  Output them.
  return I
