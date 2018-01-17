# FUNCTIONS FOR CORRELATED MUTATION ANALYSIS

import sys
import os
from math import log, fabs
from random import shuffle
from string import lower
from Bio.Clustalw import ClustalAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from numpy import *
from re import search, split


TRUE = 1
FALSE = 0



def parse_pairs_file(file_name):
  pairs = []
  f = open(file_name, 'r')
  for line in f.readlines():
    pat = search(r'^(\S+)\t(\S+)$', line)
    if pat:
      i = pat.group(1)
      j = pat.group(2)
      pairs.append((i,j))
  return pairs



def make_ordered_subalignments(aln1, aln2, pairs):
  subalign1 = ClustalAlignment()
  subalign2 = ClustalAlignment()
  for m in aln1._records:
    for n in aln2._records:
      test_pair = (m.description, n.description)
      if test_pair in pairs:
        subalign1._records.append(m)
        subalign2._records.append(n)
  return subalign1, subalign2



def condense_alignment(aln, id):
  # Condense an alignment to just those positions that are not gaps
  # in a sequence at the index given
  index = get_seq_index(aln, id)
  new_aln = ClustalAlignment()
  N = len(aln._records) # blah
  indices = range(aln.get_alignment_length())
  valid_indices = []
  # Do a filtering step to check for gaps in ref sequence 'id'
  for i in indices:
    if aln._records[index].seq.data[i] != '-':
      valid_indices.append(i)
  for n in range(N):  # For each sequences
    new_seq_str = ''
    for v in valid_indices:  # For each position where there's not a gap in 'id' sequence
      new_seq_str += aln._records[n].seq.data[v]
    new_aln._records.append(SeqRecord(Seq(new_seq_str), description=aln._records[n].description))
  return new_aln



def concat_aln_shuffle(aln1, aln2):
  concat_aln =  concatenate_alignments(aln1, aln2)
  concat_str = ''
  for k in range(len(concat_aln._records)):
    new_seq = concat_aln._records[k].seq.data 
    concat_str += new_seq + '\n'
  return concat_str



def concat_aln_sca(aln1, aln2):
  concat_aln =  concatenate_alignments(aln1, aln2)
  concat_str = ''
  for k in range(len(concat_aln._records)):
    new_seq = concat_aln._records[k].seq.data 
    desc = concat_aln._records[k].description
    concat_str += desc + '\t' + new_seq + '\n'
  return concat_str



def make_matlab_sca_script(aln_file_name, out_file_name, script_name):
  script = "[labels,aln]=get_seqs('%s');\
            [perturbation_matrix]=stat_fluc(aln);\
            [cmr]=global_sca_small(perturbation_matrix);\
            dlmwrite('%s',cmr, '\t');\
            quit;\n" % (aln_file_name, out_file_name)
  f = open(script_name, 'w')
  f.write(script)
  f.close()


  
def parse_fasta(filename):
  # Return a dictionary of sequences
  # seq_id -> sequence
  seq_dict = {}
  f = open(filename, 'r')
  for line in f.readlines():
    header_pat = search(r'^>(\S+).*$', line)
    seq_pat = search(r'^([ACDEFGHIKLMNPQRSTVWY]+)$', line)
    if header_pat:
      id = header_pat.group(1)
      seq_dict[id] = ''
    if seq_pat:
      seq_dict[id] += seq_pat.group(1)
  return seq_dict



def print_fasta(seq_dict, file_name):
  out = open(file_name, 'w')
  for s in seq_dict.keys():
    out.write(">" + s + "\n")
    out.write(seq_dict[s] + "\n")



def parse_matrix(filename):
  f = open(filename, 'r')
  data = {}
  X = []  # A list of the values on x axis
  Y = []  # list for y axis
  # Read each line
  for line in f.readlines():
    header_pat = search(r'^\t.+', line)
    data_pat = search(r'^(\S+)\t(.+)$', line)
    # If it's the first line
    if header_pat:
      array = split('\s+', line)
      for a in array:
        a_pat = search(r'^(\S+)$', a)
        if a_pat:
          X.append(a_pat.group(1))
    # Otherwise, it's data
    elif data_pat:
      y = data_pat.group(1)  # The first field tells us the current value on y axis
      Y.append(y)  # Put it in our list
      data[y] = {} # Create a data entry for y
      array = split('\t', data_pat.group(2))
      for i in range(len(array)):
        x = X[i]
        data[y][x] = float(array[i])
  # We are now finished parsing the matrix
  # Turn to a Numeric array
  I = len(X)
  J = len(Y)
  A = zeros((J,I), float)
  for x in X:
    i = X.index(x)
    for y in Y:
      j = Y.index(y)
      A[j][i] = data[y][x]
  return A, X, Y



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



def shuffle_output_to_matrix(aln1, aln2, I, J, data_file):
  N = len(I)
  data = {}
  f = open(data_file, 'r')
  for line in f.readlines():
    pat = search(r'^\((\d+),(\d+)\) : (\S+)$', line)
    if pat:
      a = int(pat.group(1))
      b = int(pat.group(2))
      p = -1 * log(float(pat.group(3)))
      if a in I and b-N in J:
        data[(a, b-N)] = p
      elif b in I and a-N in J:
        data[(b, a-N)] = p
  result = zeros((len(I), len(J)), float)
  for i in I:
    x = I.index(i)
    for j in J:
      y = J.index(j)
      if data.has_key((i,j)):
        result[x][y] = data[(i,j)]
  return result



def run_shuffle(aln1, aln2, I, J, exe, aln_file_name, out_file_name, num_cycles, num_aln):
  # Generate the output file
  cmd = "%s %s %s %d %d" % (exe, aln_file_name, out_file_name, num_cycles, num_aln)
  os.system(cmd)
  # Parse the output and put in a matrix
  N = len(I)
  data = {}
  f = open(out_file_name, 'r')
  for line in f.readlines():
    pat = search(r'^\((\d+),(\d+)\) : (\S+)$', line)
    if pat:
      a = int(pat.group(1)) - 1  # Translate to an index, 1 is really 0
      b = int(pat.group(2)) - 1
      c = float(pat.group(3))
      if c < .005:
        c = .005
      p = -1 * log(float(c))
      if a in I and b-N in J:
        data[(a, b-N)] = p
      elif b in I and a-N in J:
        data[(b, a-N)] = p
  if len(data) == 0:
    print "You have a problem: no data was found in output file"
  result = zeros((len(I), len(J)), float)
  for i in I:
    x = I.index(i)
    for j in J:
      y = J.index(j)
      if data.has_key((i,j)):
        result[x][y] = data[(i,j)]
  return result



def sca_to_matrix(aln1, aln2, I, J, data_file):
  N = len(I)
  print len(I)
  print len(J)
  result = zeros((len(I), len(J)), float)
  print result.shape
  f = open(data_file, 'r')
  # Just get the square matrix out of the file
  data = []  # Will become a list of lists
  for line in f.readlines():
    char_data = split('\t', line)
    row = []
    for c in char_data:
      row.append(float(c))
    data.append(row)
  A = array(data, float)
  # OK, now figure out which data to put in the non-square matrix
  if A.shape[0] != A.shape[1]:
    print "Your matrix is not square"
    sys.exit()
  dim = A.shape[0]
  for x in range(dim):
    for y in range(dim):
      if x in I and y-N in J:
        result[I.index(x)][J.index(y-N)] = A[x][y]
  return result



#########################
# General alignment stuff
#########################

# Shuffle the order of sequences in alignment
def shuffle_alignment(aln):
  L = len(aln._records)
  indices = range(L)
  shuffle(indices)
  new_list = []
  for i in indices:
    new_list.append(aln._records[i])
  aln._records = new_list
  return aln
 

 
####################################################################
# Functions for mapping alignment positions to residues in structure
####################################################################
def get_seq_index(aln, seq_id):  # Get the index of a particular sequence in a MSA
  print "looking for %s" % seq_id
  INDEX_FOUND = FALSE
  for a in aln._records:
    if a.description == seq_id:
      index = aln._records.index(a)
      INDEX_FOUND = TRUE
      break
  if not INDEX_FOUND:
    print "ERROR: didn't find sequence %s in list:" % seq_id
    for s in aln._records:
      print s.description,
    print
    sys.exit()
  return index



def get_aln_seq(aln, id):
  index = get_seq_index(aln, id)
  seq = aln._records[index].seq.data
  return seq



def map_aln2struc(aln, seq):
  # This is a ridiculously simple function. All it does is give you a mapping
  # of sequence residue indices to structure residue indices. Since there are
  # possibly gaps in the sequence but never in the structure (if it's a single
  # polypeptide, which it should be), all you have to do is not increment the
  # structure index if the sequence has a gap

  # This function assumes you know that the sequence from the alignment, passed
  # as argument 'seq', is identical to that in the pdb file
  aln2struc = {}
  x = 0
  for i in range(len(seq)):
    if seq[i] != '-':  # If the residue is anything but a gap
      aln2struc[i] = x
      x += 1
    else:
      aln2struc[i] = None  # May want to rethink this at some point
  return aln2struc



# Take two seq_dicts and a list of pairs of seq_ids [(m0, n0), (m1, n1), ...] and return two edited
# seq_dicts such that one list of keys is [m0, m1, ...], the other is [n0, n1, ...]
def edit_sequence_dicts(seq_dict1, seq_dict2, pairs):
  # Edit the first seq_dict
  for m in seq_dict1.keys():
    MATCH_FOUND = FALSE
    for i,j in pairs:
      # Figure out if this m is equal to some i in the pairs
      if m == i:
        MATCH_FOUND = TRUE
    if not MATCH_FOUND:
      del seq_dict1[m]
  # Do it again for the other seq_dict
  for n in seq_dict2.keys():
    MATCH_FOUND = FALSE
    for i,j in pairs:
      # Figure out if this m is equal to some i in the pairs
      if n == j:
        MATCH_FOUND = TRUE
    if not MATCH_FOUND:
      del seq_dict2[n]



###############################################
# Analyze alignment within context of structure
###############################################
def valencia_fly(aln1, aln2, I, J, subs_mat, COMPLEX_FLAG):
  N = len(aln1._records) 
  result = zeros((len(I), len(J)), float)

  # GET THE MEANS AND SDS OF ALL COLUMNS FIRST
  # THESE ARE COSTLY TO COMPUTE
  mean_i = {}
  mean_j = {}
  sd_i = {}
  sd_j = {}

  # MEANS
  sum = 0.0
  for i in I:
    m = aln1.get_column(i)
    for l in m:
      if l != '-':
        for k in m:
	  if k != '-':
	    sum += subs_mat[l][k]
    mean_i[i] = sum / N
  
  sum = 0.0
  for j in J:
    n = aln2.get_column(j)
    for l in n:
      if l != '-':
	for k in n:
	  if k != '-':
	    sum += subs_mat[l][k]
    mean_j[j] = sum / N

  # SD
  sum = 0.0
  for i in I:
    m = aln1.get_column(i)
    for l in m:
      for k in m:
        if l != '-' and k != '-':
          sum += pow(subs_mat[l][k] - mean_i[i], 2)
        else:
          sum += pow(0.0 - mean_i[i], mean_i[i], 2)
    sd_i[i] = pow(sum / N, .5)
  
  sum = 0.0
  for j in J:
    n = aln2.get_column(j)
    for l in n:
      for k in n:
        if l != '-' and k != '-':
          sum += pow(subs_mat[l][k] - mean_j[j], 2)
        else:
          sum += pow(0.0 - mean_i[i], mean_j[j], 2)
    sd_j[j] = pow(sum / N, .5)

  # FOR EACH PAIR OF COLUMNS 
  for i in I: 
    m = aln1.get_column(i)
    x = I.index(i)
    for j in J:  
      n = aln2.get_column(j)
      y = J.index(j)
      sum = 0.0
            

def valencia(aln1, aln2, I, J, subs_mat, w, COMPLEX_FLAG):
  n_off = 0  # The number of positions off the diagonal that a pair of residues must be.
  if not COMPLEX_FLAG:
    print "Excluding pairs less than %d positions apart in sequence" % n_off
             # This prevents finding correlations between residues nearby in sequence
	     # and thus biasing the results.
	     # Set to 0 to analyze all pairs including self-self
	     # Set to 1 to exclude self-self
	     # Set to 2 to exclude those 1 off the diagonal
	     # Etc.
  N = len(aln1._records) 
  print "Converting alignment to matrix format"
  sys.stdout.flush()
  S1, MU1, SIGMA1 = get_valencia_matrices(I, aln1, subs_mat)
  S2, MU2, SIGMA2 = get_valencia_matrices(J, aln2, subs_mat)
  result = zeros((len(I), len(J)), float)
  print "Analyzing positional correlations"
  sys.stdout.flush()
  for i in I:  # for each columm
    x = I.index(i)
    for j in J:  # for each column
      y = J.index(j)
      if abs(x-y) < n_off and not COMPLEX_FLAG:
        continue  # skip to the next i,j
      else:
	sys.stdout.flush()
        sum = 0.0  # Sum in numerator in correlation
        for k in range(N):  # for each amino acid
          for l in range(N):  # for each amino acid
            # Just do correlation here..., don't call a function
	    # DEBUG
	    a = w[k][l]
	    b = S1[x][k][l]
	    c = MU1[x]
	    d = S2[y][k][l]
	    e = MU2[y]
	    sum +=  (a * (b - c) * (d - e)) 
        result[x][y] = abs(sum / (N * N * SIGMA1[x] * SIGMA2[y]))
  return result



def valencia_dimer_permute(aln1, aln2, I, J, subs_mat, w):
  num_permutations = 100
  num_sig = 0
  num_ins = 0
  N = len(aln1._records) 
  print "Permuting for %d iterations" % num_permutations
  S1, MU1, SIGMA1 = get_valencia_matrices(I, aln1, subs_mat)
  S2, MU2, SIGMA2 = get_valencia_matrices(J, aln2, subs_mat)
  result = zeros((len(I), len(J)), float)
  for i in I:  # for each columm
    x = I.index(i)
    for j in J:  # for each column
      y = J.index(j)
      print i, j
      sum = 0.0  # Sum in numerator in correlation
      for k in range(N):  # for each amino acid
        for l in range(N):  # for each amino acid
          # Just do correlation here..., don't call a function
	  # DEBUG
	  a = w[k][l]
	  b = S1[x][k][l]
	  c = MU1[x]
	  d = S2[y][k][l]
	  e = MU2[y]
	  sum +=  (a * (b - c) * (d - e)) 
      score = sum / (N * N * SIGMA1[x] * SIGMA2[y])
      # Now permute rows many times and see if score is still this high
      permuted_scores = zeros(num_permutations, float)  # An empty vector
      m = aln1.get_column(i)
      for p in range(num_permutations):
        n = randomize_string(aln2.get_column(j))
	q = valencia_correlation(m,n,subs_mat,w) 
        permuted_scores[p] = q
      lower, upper = r.quantile(permuted_scores, [0.05, 0.95])
      if score <= lower or score >= upper:
        result[x][y] = abs(score)  # Here is where we do the Fodor and Aldrich abs(r) thing
	num_sig += 1
      else:
        result[x][y] = 0.0
	num_ins += 1
  print "%f of scores were significant" % (num_sig / float(num_sig + num_ins))
  print "%f of scores were not" % (num_ins / float(num_sig + num_ins))
  return result
 

 
def get_valencia_indices(aln, entropy_threshold, gap_threshold):
  indices = []
  L = aln.get_alignment_length()
  for i in range(L):
    col = aln.get_column(i)
    ################################################
    # Check that the column is totally conserved
    ################################################
    h = entropy(col)
    if h <= entropy_threshold:  # If the column is totally conseved
      continue
    #################################################
    # Check that the column is more than 10% gaps
    #################################################
    gap_count = 0
    N = len(col)
    for a in col:
      if a == '-':
        gap_count += 1
    if gap_count / float(N) > gap_threshold:
      continue
    indices.append(i)
  return indices


  
def get_valencia_structure_indices(aln, index, entropy_threshold, gap_threshold):
  indices = []
  L = aln.get_alignment_length()
  for i in range(L):
    m = aln.get_column(i)
    if column_is_analyzable(m, index, entropy_threshold, gap_threshold):
      indices.append(i)
    else:
      pass
  return indices



def column_is_analyzable(col, index, entropy_threshold, gap_threshold):
  ################################################
  # If column is a gap in structure sequence, 
  # forget about it
  ################################################
  if col[index] == '-':
    #print "Deleted column because structure contains a gap"
    #print col
    return FALSE
  ################################################
  # Check that the column is totally conserved
  ################################################
  h = entropy(col)
  if h <= entropy_threshold:  # If the column is totally conseved
    return FALSE
  #################################################
  # Check that the column is more than 10% gaps
  #################################################
  gap_count = 0
  N = len(col)
  for a in col:
    if a == '-':
      gap_count += 1
  if gap_count / float(N) > gap_threshold:
    return FALSE
  # If we get to this point the column is OK
  return TRUE



def get_omes_indices(aln, index, gap_threshold):
  indices = []
  L = aln.get_alignment_length()
  for i in range(L):
    m = aln.get_column(i)
    if column_is_analyzable(m, index, 0.0, gap_threshold):
      indices.append(i)
  return indices



def get_mi_indices(aln, index):
  indices = []
  L = aln.get_alignment_length()
  for i in range(L):
    m = aln.get_column(i)
    if m[index] != '-':  # If there's not a gap in the structure 
      indices.append(i)
  return indices
 

 
def get_sca_indices(aln, index):
  indices = []
  L = aln.get_alignment_length()
  for i in range(L):
    m = aln.get_column(i)
    if m[index] != '-':  # If there's not a gap in the structure 
      indices.append(i)
  return indices



def get_indices(aln, index):
  indices = []
  L = aln.get_alignment_length()
  for i in range(L):
    m = aln.get_column(i)
    if m[index] != '-':  # If there's not a gap in the structure 
      indices.append(i)
  return indices


  
def get_indices2(aln, index):
  L = aln.get_alignment_length()
  ngap = 0
  for i in range(L):
    m = aln.get_column(i)
    if m[index] == '-':  # If there's not a gap in the structure 
      ngap += 1
  return range(L-ngap)


  
def omes(aln1, aln2, I, J, COMPLEX_FLAG):
  # Method of Kaas and Horovitz, Proteins 2002
  # Observed Minus Expected Squared covariance algorithm
  # A data structure for the result
  result = zeros((len(I), len(J)), float)
  for i in I:
    for j in J:
      if i == j and not COMPLEX_FLAG:
          continue  # skip to the next i,j
      else:
        sys.stdout.flush()
        a = I.index(i)
        b = J.index(j)
        m = aln1.get_column(i)
        n = aln2.get_column(j)
        ################################
        # Calculate the score here
        #################################
        pair_dict = {}   # Dictionary of all observed pairs of amino acids
        aa_dict_i = {}   # Dictionary of counts of kinds of amino acids in column i
        aa_dict_j = {}   # Same for column j
        Nvalid = 0
        # What amino acid pairs are observed?
        for k in range(len(m)):
          # Extract the amino acids from columns
          x = m[k]
	  y = n[k]
	  # Deal with pairs
	  if x != '-' and y != '-':  # It is a valid pair
	    Nvalid += 1
	    if not pair_dict.has_key((x,y)):
	      pair_dict[(x,y)] = 1   # We have observed the first pair. Values in dict are Nex.
	    else:
	      pair_dict[(x,y)] += 1  # Add this observed count
	  # Deal with individual amino acid occurrences
          if x != '-':
	    if not aa_dict_i.has_key(x):
	      aa_dict_i[x] = 1
	    else:
	      aa_dict_i[x] += 1
          if y != '-':
	    if not aa_dict_j.has_key(y):
	      aa_dict_j[y] = 1
	    else:
	      aa_dict_j[y] += 1
        # Calculate the terms of the score for each observed amino acid pair
        score = 0.0  # This will be a sum of scores for each amino acid pair
        for x, y in pair_dict.keys():
          Cxi = aa_dict_i[x]
	  Cyj = aa_dict_j[y]
  	  Nex = Cxi * Cyj / float(Nvalid)
	  Nobs = pair_dict[(x,y)]
	  #  CHOOSE AN EQUATION:
	  # This is the equation reported in Fodor's paper
	  #score += pow(Nobs - Nex, 2) / float(Nvalid)
	  # This is the version from the original paper, Kaas & Horovitz
	  score += pow(Nobs - Nex, 2) / float(Nex)
	  # score is the chi squared statistic
        df = (len(aa_dict_i) * len(aa_dict_j)) - 1
        pval = 1-r.pchisq(score, df)
	
	###########################################
	# Only store significant scores given df:
        #if pval <= .05: 
        #  result[a][b] = score  # Store the result
	###########################################
	
        result[a][b] = score  # Store the result
  return result



def mi_old(aln1, aln2, X, Y, COMPLEX_FLAG):
  # The mutual information method of Atchley et al.
  result = zeros((len(X), len(Y)), float)
  for x in X:
    for y in Y:
      if x == y and not COMPLEX_FLAG:
	  continue  # Skip analyzing the self-self mutual information
      else: 
	sys.stdout.flush()
        a = X.index(x)
        b = Y.index(y)
        col1 = aln1.get_column(x)
        col2 = aln2.get_column(y)
        pair_dict = {}
        aa_dict_x = {}
        aa_dict_y = {}
        Nvalid = 0
        score = 0.0
        for i in range(len(col1)):
          m = col1[i]
	  n = col2[i]
	  if m != '-' and n != '-':
	    Nvalid += 1
	    if not pair_dict.has_key((m,n)):
	      pair_dict[(m,n)] = 1
	    else:
	      pair_dict[(m,n)] += 1
            if not aa_dict_x.has_key(m):
	      aa_dict_x[m] = 1
	    else:
	      aa_dict_x[m] += 1
            if not aa_dict_y.has_key(n):
	      aa_dict_y[n] = 1
	    else:
	      aa_dict_y[n] += 1
        for j,k in pair_dict.keys():
          pjk = pair_dict[(j,k)] / float(Nvalid)
  	  pj = aa_dict_x[j] / float(Nvalid)
  	  qk = aa_dict_y[k] / float(Nvalid)
          if pjk > 0.0:
	    score += pjk * log(pjk / (pj * qk))
	
        result[a][b] = score  # Store the result
  return result



def mi(aln1, aln2, I, J, COMPLEX_FLAG, NORM_FLAG):
  # The normalized mutual information method of LM Wahl, Bioinformatics 2005
  result = zeros((len(I), len(J)), float)
  for i in I:
    for j in J:
      if i == j and not COMPLEX_FLAG:
	  continue  # Skip analyzing the self-self mutual information
      else: 
        a = I.index(i)
        b = J.index(j)
        m = aln1.get_column(i)
        n = aln2.get_column(j)
        score = mutual_information(m,n)
	if NORM_FLAG:
	  if joint_entropy(m,n) > 0.0:
	    score /= joint_entropy(m,n)
        result[a][b] = score  # Store the result
    sys.stdout.flush()
  return result
  


def normalized_mi(aln1, aln2, I, J, COMPLEX_FLAG):
  # The normalized mutual information method of LM Wahl, Bioinformatics 2005
  result = zeros((len(I), len(J)), float)
  for i in I:
    for j in J:
      if i == j and not COMPLEX_FLAG:
	  continue  # Skip analyzing the self-self mutual information
      else: 
        a = I.index(i)
        b = J.index(j)
        m = aln1.get_column(i)
        n = aln2.get_column(j)
	if joint_entropy(m,n) > 0.0:
          score = mutual_information(m,n) / joint_entropy(m,n)
          result[a][b] = score  # Store the result
  return result


  
def entropy_dimer(aln1, aln2, I, J):
  # A data structure for the result
  result = zeros((len(I), len(J)), float)
  
  for i in I:
    for j in J:
      x = I.index(i)
      y = J.index(j)
      m = aln1.get_column(i)
      n = aln2.get_column(j)
      H_m = entropy(m)
      H_n = entropy(n)
      r = entropy(m + n)
      result[x][y] = r  # Store the result
  return result



def normalize_result(M):
  max_score = 0.0
  for row in M:
    for item in row:
      if item > max_score:
        max_score = item
  print "max score: %f" % max_score
  M /= max_score



############################################
# Functions for determining residue distance
# within or between proteins
############################################
def get_distance(A, B, struc, I, J):
  # We aren't thinking about gaps in an alignment or how they relate to structure
  distance = zeros((len(I), len(J)), float)
  for i in I:
    for j in J:
      a = struc[0][A].child_list[i]  # Should result in a residue	
      b = struc[0][B].child_list[j]  # Should result in a residue	
      d = get_residue_distance(a, b)  # Distance between residues x and y
      distance[i][j] = d
  return distance



def get_distance_subset(A, B, struc, aln_to_struc, I, J):
  # I and J are lists of indices of alignment columns that are valid to
  # analyze by some criterion.
  # Make a data structure for the result
  distance = zeros((len(I), len(J)), float)
  for i in I:
    for j in J:
      # Translate the indices
      x = I.index(i)
      y = J.index(j)
      # Get the distance
      a = struc[0][A].child_list[aln_to_struc[A][i]]  # Should result in a residue	
      b = struc[0][B].child_list[aln_to_struc[B][j]]  # Should result in a residue	
      d = get_residue_distance(a, b)  # Distance between residues x and y
      distance[x][y] = d
  return distance



def get_distance_subset_monomer(chain_id, struc, aln_to_struc, I, J):
  # I and J are lists of indices of alignment columns that are valid to
  # analyze by some criterion.
  # Make a data structure for the result
  distance = zeros((len(I), len(J)), float)
  for i in I:
    for j in J:
      # Translate the indices
      x = I.index(i)
      y = J.index(j)
      # Get the distance
      a = struc[0][chain_id].child_list[aln_to_struc[chain_id][i]]  # Should result in a residue	
      b = struc[0][chain_id].child_list[aln_to_struc[chain_id][j]]  # Should result in a residue	
      d = get_residue_distance(a, b)  # Distance between residues x and y
      distance[x][y] = d
  return distance



def get_distance_dimer(a1, a2, A, B, struc, aln_to_struc):
  # What are the distances between all residues in the structure?
  # Determine alignment lengths
  # Exclude gap positions
  n_gaps_i = 0
  n_gaps_j = 0
  for a in a1:
    if a == '-':
      n_gaps_i += 1
  for a in a2:
    if a == '-':
      n_gaps_j += 1
  I = len(a1)
  J = len(a2)
  # Make a data structure for the result
  distance = zeros((I,J), float)
  for i in range(I):
    for j in range(J):
      if a1[i] != '-' and a2[j] != '-':  # If neither position is a gap    
	# Get the distance
        x = struc[0][A].child_list[aln_to_struc[A][i]]  # Should result in a residue	
        y = struc[0][B].child_list[aln_to_struc[B][j]]  # Should result in a residue	
	d = get_residue_distance(x, y)  # Distance between residues x and y
        distance[i][j] = d
  return distance



def get_distance_monomer(struc):
  # What are the distances between all residues in the structure?
  # Determine alignment lengths
  # Exclude gap positions
  # Make a data structure for the result
  L = len(struc[0].child_list[0].child_list)
  distance = zeros((L,L), float)
  for i in range(L):
    for j in range(L):
      # Get the distance
      x = struc[0].child_list[0].child_list[i]
      y = struc[0].child_list[0].child_list[j]
      d = get_residue_distance(x, y)  # Distance between residues x and y
      distance[i][j] = d
  return distance



def get_residue_distance(x, y):
  # USE CB IF POSSIBLE
  if x.child_dict.has_key('CB'):
    atom1 = x['CB']
  else:
    atom1 = x['CA']
  if y.child_dict.has_key('CB'):
    atom2 = y['CB']
  else:
    atom2 = y['CA']
  return atom1 - atom2



def get_closest_distance(x,y):
  # Get the shortest distance between any two atoms of a pair of residues
  min = x['CA'] - y['CA']  # We need an initial value for min
  for i in x.child_list:
    for j in y.child_list:
      d = i - j
      if d < min:
        min = d
  return min


	
def filter_result_by_distance(dist, corr, max_d):
  if len(dist) != len(corr):
    print "Vectors different lengths"
    sys.exit()
  N = len(dist)
  X = []
  Y = []
  for i in range(N):
    if dist[i] <= max_d:
      X.append(dist[i])
      Y.append(corr[i])
  return X, Y



def charged_pairs(col1, col2):
  plus = ['K', 'R']
  minus = ['D', 'E']
  
  if len(col1) != len(col2):
    print "Unequal column lengths"
    sys.exit()
  N = len(col1)
  plus_minus = 0
  minus_plus = 0
  for i in range(N):
    a = col1[i]
    b = col2[i]
    if a in plus and b in minus:
      plus_minus += 1
    elif b in plus and a in minus:
      minus_plus += 1
  return plus_minus * minus_plus


    
#############################
# Output formatting functions
#############################
def write_matrix(A, X, Y, filename):
  f = open(filename, 'w')
  # Print out the top row
  for j in range(len(Y)):
    f.write("\t%s" % Y[j])
  f.write("\n")
  # Print out the data
  for i in range(len(X)):
    f.write("%s" % X[i])
    for j in range(len(Y)):
      f.write("\t%f" % A[i][j])
    f.write("\n")
  f.close()
    


def write_dict_matrix(I, J, mat, outfile):
  # print a header line
  for i in I:
    outfile.write("\t%s" % i)
  outfile.write("\n")
  for j in J:
    outfile.write("%s" % j)
    for i in I:
      outfile.write("\t")
      outfile.write("%f" % mat[i][j])
    outfile.write("\n")
  


def write_result_matrix(I, J, mat, outfile):
  # print a header line
  for i in I:
    outfile.write("\t%d" % i)
  outfile.write("\n")
  for j in J:
    y = J.index(j)
    outfile.write("%d" % j)
    for i in I:
      x = I.index(i)
      outfile.write("\t")
      outfile.write("%f" % mat[x][y])
    outfile.write("\n")



def write_protein_result_matrix(I, J, s, A, B, mat, outfile):
  # print a header line
  for i in I:
    res_i = s[0][A].child_list[i].get_resname() + str(s[0][A].child_list[i].id[1])
    outfile.write("\t(%d)%s" % (i, res_i))
  outfile.write("\n")
  for j in J:
    res_j = s[0][B].child_list[j].get_resname() + str(s[0][B].child_list[j].id[1])
    outfile.write("(%d)%s" % (j,res_j))
    for i in I:
      outfile.write("\t")
      outfile.write("%f" % mat[i][j])
    outfile.write("\n")



def write_column_result_list(aln1, aln2, s, A, B, distance, result, outfile):
  for i in range(aln1.get_alignment_length()):
    res_i = s[0][A].child_list[i].get_resname() + str(s[0][A].child_list[i].id[1])
    col_i = aln1.get_column(i)
    for j in range(aln2.get_alignment_length()):
      res_j = s[0][B].child_list[j].get_resname() + str(s[0][B].child_list[j].id[1])
      col_j = aln2.get_column(j)
      outfile.write("%f\t%f\t%d\t%s\t%s\tx\t%d\t%s\t%s\n" % (result[i][j], distance[i][j], i, res_i, col_i, j, res_j, col_j))



def write_protein_result_list(I, J, s, A, B, distance, result, outfile):
  # print a header line
  for j in J:
    res_j = s[0][B].child_list[j].get_resname() + str(s[0][B].child_list[j].id[1])
    for i in I:
      res_i = s[0][A].child_list[i].get_resname() + str(s[0][A].child_list[i].id[1])
      outfile.write("%s\t%s\t%f\t%f\n" % (res_i, res_j, distance[i][j], result[i][j]))



def write_result_table(A, B, outfile):
  # Do the matrices have appropriate dimensions?
  # We want them to be the same size
  if A.shape != B.shape:
    print "ERROR: matrices of different dimensions"
    sys.exit()
  if len(A.shape) > 1:
    print "ERROR: your output matrices are not unit vectors"
    sys.exit()
  N = len(A)
  for i in range(N):
    outfile.write("%f\t%f\n" % (A[i], B[i]))



def write_correlation_to_structure(result, struc, aln_to_struc, A, B, aln1, aln2, I, J, index1, index2):
  # For each pair of columns in our alignments, what are the residues
  # in the known structure that correspond to these columns. So if we
  # see a highly correlated pair of positions, we know which residues 
  # to look at.
  # Output something like:
  #   r res1 res2
  # or:
  #   r res1 res2 col1 col2
  for i in I:
    for j in J:
      x = I.index(i)
      y = J.index(j)
      r = result[x][y]
      a = struc[0][A].child_list[aln_to_struc[A][i]]  # Should result in a residue	
      b = struc[0][B].child_list[aln_to_struc[B][j]]  # Should result in a residue	
      m = convert_column_case(aln1.get_column(i), index1)
      n = convert_column_case(aln2.get_column(j), index2)
      d = get_residue_distance(a, b)  # Distance between residues x and y
      print "%s\t%s\t%s\t%s\t%f\t%f" % (a.resname + str(a.id[1]), m, b.resname + str(b.id[1]), n, r, d)



def convert_column_case(column, index):
  a = column[:index]
  b = lower(column[index])
  c = column[index + 1:]
  return a + b + c
  


##################################################
# Statistics and other numerical calculations
##################################################
def entropy(X):
  base = 20
  alphabet = {}
  for x in X:
    if x not in alphabet.keys():
      alphabet[x] = 1
    else:
      alphabet[x] += 1
  L = len(X)
  N = len(alphabet.keys())
  H = 0.0
  for a in alphabet.keys():
    p = alphabet[a] / float(L)
    H -= p * log(p, base)
  return H



def joint_entropy(X, Y):
  if len(X) != len(Y):
    print 'ERROR: vectors different lengths'
    sys.exit()
  base = 20
  alphabet = {}
  L = len(X)
  for i in range(L):
    x = X[i]
    y = Y[i]
    if (x,y) not in alphabet.keys():
      alphabet[(x,y)] = 1
    else:
      alphabet[(x,y)] += 1
  N = len(alphabet.keys())
  H = 0.0
  for a in alphabet.keys():
    p = alphabet[a] / float(L)
    H -= p * log(p, base)
  return H



def mutual_information(X,Y):
  Hx = entropy(X)
  Hy = entropy(Y)
  Hxy = joint_entropy(X,Y)
  MI = Hx + Hy - Hxy
  return MI



def fishers_test_for_binary_matrix_overlap(A, B):
  m = 0
  n = 0
  k = 0
  N = 0
  for i in range(len(A)):
    for j in range(len(A[i])):
      N += 1
      if A[i][j] and B[i][j]:
        k += 1
      elif A[i][j]:
        m += 1
      elif B[i][j]:
        n += 1
  return k, m, n, N



def fishers_exact_test(A,B):
  # Takes two binary vectors as args
  # Get m, n, k, and N
  m = 0
  n = 0
  k = 0
  N = 0
  # Check that we have same length vectors
  if len(A) != len(B):
    print "ERROR: Unequal length vectors"
    return
  for i in range(len(A)):
      N += 1
      if A[i] and B[i]:
        k += 1
        m += 1
        n += 1
      elif A[i]:
        m += 1
      elif B[i]:
        n += 1
  print m, n, k, N
  # A function for factorials
  f = fact
  pval =  calc_hyper_pval(k, m, n, N, fact)
  return pval



def calc_hyper_p(k, m, n, N, fact):
  ################################################################
  # Get the probability density from a hypergeometric distribution
  # at a particular k
  ################################################################
  
  #numerator = fact(m) * fact(n) * fact(N-m) * fact(N-n)
  #denominator = fact(k) * fact(n-k) * fact(m-k) * fact(N-m-n+k)
  #p = float(numerator) / float(denominator)
  p = r.choose(m,k) * r.choose(N-m, n-k) / float(r.choose(N,n))
  return p



def calc_hyper_pval(k, m, n, N, fact):
  ################################################################
  # Get a P-value, summing over a range of possible k, using
  # the probability density from a hypergeometric distribution
  # at each particular k
  ################################################################
  sum = 0.0
  for x in range(k, min(m,n)+1):
    p = calc_hyper_p(x, m, n, N, fact)
    print x, p
    sum += p
  return sum



def fact(n):
  # Calculate n!
  if n > 170:  # The max on this (or most) machines
    a = pow(2 * 3.14159 * n, .5) * pow(n / exp(1), n)
  else:
    a =  r.gamma(n+1)
  return a


# Randomize a matrix.
def randomize_vector(A):
  B = zeros(len(A), float)
  list = []
  for a in A:
    list.append(a)
  shuffle(list)
  for i in range(len(list)):
    B[i] = list[i]
  return B



def randomize_sequence(seq):  # Takes a Bio.SeqRecord.SeqRecord object
  res_list = []
  for a in seq.seq.data:
    res_list.append(a)
  shuffle(res_list)
  new_string = ''
  for a in res_list:
    new_string += a
  seq.seq.data = new_string


  
def randomize_string(str):
  char_list = []
  for a in str:
    char_list.append(a)
  shuffle(char_list)
  new_string = ''
  for a in char_list:
    new_string += a
  return new_string


  
def get_subs_mat(file_name):
  ###########################################################
  # Read an amino acid substitution matrix from at text file
  # in the usual format.
  # Return a dictionary of dictionaries corresponding to rows
  # and columns in the matrix. 
  # Crash if the rows and columns indices are not identical
  ###########################################################
  f = open(file_name, 'r')
  subs_mat = {}
  KEYS_READ = FALSE
  for line in f.readlines():
    if line[0] != '#':
      if not KEYS_READ:
        fields = split('\s+', line)
	KEYS_READ = TRUE
	keys = []
	for k in fields:
	  if k != '':
	    keys.append(k)
	    subs_mat[k] = {}
      elif KEYS_READ:
        array = split('\s+', line)
	k = array[0]
	values = array[1:-1]
	if len(values) != len(keys):
	  print "We got a problem: keys and values have different lengths: %d %d" % (len(keys), len(values))
	  sys.exit()
	else:
	  for i in range(len(keys)):
	    subs_mat[k][keys[i]] = float(values[i])
  return subs_mat



def print_subs_mat(mat):
  print ' ',
  for i in mat.keys():
    print "%2s" % i,
  print
  for i in mat.keys():
    print "%s" % i,
    for j in mat[i]:
      print "%2s" % mat[i][j],
    print



########################################################################
# Functions for analyzing correlated mutations as in the method of Gobel
########################################################################
def valencia_correlation(col1, col2, subs_mat, w):
  # Check that columns are same length
  if len(col1) != len(col2):
    print "Columns unequal length"
    sys.exit()
  N = len(col1)
  # Make the position-specific substitution matrices
  mat1 = zeros((N,N), float)
  mat2 = zeros((N,N), float)
  for k in range(N):
    for l in range(N):
      m = col1[k]
      n = col1[l]
      if m == '-' or n == '-':
        mat1[k][l] = 0.0
      else:
        mat1[k][l] = subs_mat[m][n]
  for k in range(N):
    for l in range(N):
      m = col2[k]
      n = col2[l]
      if m == '-' or n == '-':
        mat2[k][l] = 0.0
      else:
        mat2[k][l] = subs_mat[m][n]
  r = calc_weighted_matrix_correlation(mat1, mat2, w)
  return r  # Fodor and aldrich did this



def valencia_shuffle(aln1, aln2, I, J, subs_mat):
  result = zeros((len(I), len(J)), float)
  for i in I:
    col1 = aln1.get_column(i)
    for j in J:
      col2 = aln2.get_column(j)
      print i, '\t', col1
      print j, '\t', col2
      p = valencia_pval_clever(col1, col2, subs_mat)
      result[i][j] = p  
      print p
      print
  return result



def valencia_pval(col1, col2, subs_mat):
  n_trials = 1000
  # Check that columns are same length
  if len(col1) != len(col2):
    print "Columns unequal length"
    sys.exit()
  # DO THE REFERENCE CALCULATION
  # Make the position-specific substitution matrices
  mat1 = valencia_pssm(col1, subs_mat)
  mat2 = valencia_pssm(col2, subs_mat)
  r = abs(calc_matrix_correlation(mat1, mat2))
  n_greater = 0
  for trial in range(n_trials):
    new_col = randomize_string(col2)
    new_mat = valencia_pssm(new_col, subs_mat)
    test_r = abs(calc_matrix_correlation(mat1, new_mat))
    if test_r >= r:
      n_greater += 1
  pval = float(n_greater) / n_trials
  return pval



def valencia_pval_clever(col1, col2, subs_mat):
  if entropy(col1) == 0.0 or entropy(col2) == 0.0:
    pval = 1.0
  else:
    n_trials = (10, 100, 1000)
    # Check that columns are same length
    if len(col1) != len(col2):
      print "Columns unequal length"
      sys.exit()
    # DO THE REFERENCE CALCULATION
    # Make the position-specific substitution matrices
    mat1 = valencia_pssm(col1, subs_mat)
    mat2 = valencia_pssm(col2, subs_mat)
    r = abs(pearsonr(mat1.flatten(), mat2.flatten())[0])
    for n in n_trials:
      n_asgood = 0
      for trial in range(n):
        new_col = randomize_string(col2)
        new_mat = valencia_pssm(new_col, subs_mat)
        test_r = abs(pearsonr(mat1.flatten(), new_mat.flatten())[0])
        if test_r >= r:
          n_asgood += 1
      if n_asgood > 0:
        pval = float(n_asgood) / n
        break
    if n_asgood == 0:
      pval = 0.0
      pval = 1 / float(n_trials[-1])   # e.g., 1/1000
  return pval



def valencia_pssm(col, subs_mat):
  N = len(col)
  mat = numpy.zeros((N,N), numpy.float)
  for k in range(N):
    for l in range(N):
      m = col[k]
      n = col[l]
      if m == '-' or n == '-':
        mat[k][l] = 0.0
      else:
        mat[k][l] = subs_mat[m][n]
  return mat


  
def get_valencia_matrices(I, aln, subs_mat):
  # I is the list of alignment indices that are valid to analyze
  # Turn a column from an alignment of N sequences into an N x N substitution matrix
  # Also return the mean and sd for each matrix
  # So return a list of 3 things for each column: matrices, mean, sd
  # 3 lists of length L
  
  N = len(aln.get_column(0))  # Figure out how many sequences there are
  L = len(I)  # This time L is just the number of columns that can be analyzed. 
              # So it is potentially less than the full alignment length
  
  # Define the data structures returned
  S = zeros((L,N,N), float)  # All the PSSMs for each position
  MU = zeros(L, float)
  SIGMA = zeros(L, float)
  
  for i in I:
    x = I.index(i)
    col = aln.get_column(i)
    # Make the position-specific substitution matrices
    for k in range(N):
      for l in range(N):
        m = col[k]
        n = col[l]
        #if m == '-' or n == '-':
        if m == '-' or n == '-' or m == 'X' or n == 'X':
          S[x][k][l] = 0.0
        else:
          S[x][k][l] = subs_mat[m][n]
    # Get the mean and sd
    # To do this make a temporary unit vector of values
    v = matrix_to_unit_vector(S[x])
    MU[x] = mean(v)
    SIGMA[x] = sd(v)
  return S, MU, SIGMA

  

###################################################################
# Functions for analyzing statistical coupling of residue positions
# in a MSA as in Lockless and Ranganathan, Science 1999
###################################################################

def sca_monomer(aln, aa_dict, struc, aln_to_struc, chain_id):
  kT = 1.0  # An arbitrary energy unit. Could think about changing this.
  # List of amino acids to consider
  aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
  # The length of the alignment
  L = aln.get_alignment_length()
  # The number of sequences
  N = len(aln._records)
  
  # Get vectors of amino acid counts for each column
  F = get_aa_counts(L, aln, aa_list)   # A list of 20-long unit vectors of aa counts for each alignment column

  # Get a vector of amino acid frequencies for all proteins everywhere
  # This just comes from supplied aa_dict
  A = zeros(20,float)
  for a in aa_dict.keys():
    if a in aa_list:
      x = aa_list.index(a)
      A[x] = aa_dict[a]
  # DEBUG

  # Now get the binomial densities P(x) for each column, each amino acid
  P = get_P(L, aln, aa_list, F, N, A) 

  # Get the counts vector for a hypothetical column in which all amino acids are observed the 
  # expected number of times given their frequency in the MSA
  H = get_H(L, aln, aa_list)
  
  # Get P_MSA
  # Vector of amino acid frequencies in MSA
  P_MSA = get_P_MSA(N, A, H)

  
  # Get vectors of statistical energies
  # For each column, for each amino acid
  G_stat_ix = get_G_stat_ix(L, P, aa_list, P_MSA, kT)

  # Get the vectors of statistical energies
  # for each alignment position, summed over all 20 amino acids
  # A measure of how out of equilibrium each alignment position is
  G_stat_i = get_G_stat_i(L, P, P_MSA, kT)
  
  for i in aln_to_struc[chain_id].keys():  # For a bunch of integers representing alignment column indices
    s = aln_to_struc[chain_id][i]          # get another integer, representing the index of the residue in structure
    if s:
      t = struc[0][chain_id].child_list[s]   # get the residue object for that index
      res_name = t.resname + str(t.id[1])
      #print "%s\t%f" % (res_name, G_stat_i[i])
  result = perturb(aln, aa_list, P, P_MSA, L, N, A, kT)
  #print result
  return result


  
def sca_dimer(aln1, aln2, aa_dict, s, aln_to_struc, I, J, chain_id1, chain_id2):
  kT = 1.0  # An arbitrary energy unit. Could think about changing this.
  # List of amino acids to consider
  aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
  # The length of the alignment
  L = aln1.get_alignment_length()
  # The number of sequences
  N = len(aln1._records)
  
  # Get vectors of amino acid counts for each column
  F1 = get_aa_counts(L, aln1, aa_list)   # A list of 20-long unit vectors of aa counts for each alignment column
  F2 = get_aa_counts(L, aln2, aa_list)   # A list of 20-long unit vectors of aa counts for each alignment column

  # Get a vector of amino acid frequencies for all proteins everywhere
  # This just comes from supplied aa_dict
  A = zeros(20,float)
  for a in aa_dict.keys():
    if a in aa_list:
      x = aa_list.index(a)
      A[x] = aa_dict[a]
  # DEBUG

  # Now get the binomial densities P(x) for each column, each amino acid
  P1 = get_P(L, aln1, aa_list, F1, N, A) 
  P2 = get_P(L, aln2, aa_list, F2, N, A) 

  # Get the counts vector for a hypothetical column in which all amino acids are observed the 
  # expected number of times given their frequency in the MSA
  H1 = get_H(L, aln1, aa_list)
  H2 = get_H(L, aln2, aa_list)
  
  # Get P_MSA
  # Vector of amino acid frequencies in MSA
  P_MSA1 = get_P_MSA(N, A, H1)
  P_MSA2 = get_P_MSA(N, A, H2)
  
  # Get vectors of statistical energies
  # For each column, for each amino acid
  G_stat_ix1 = get_G_stat_ix(L, P1, aa_list, P_MSA1, kT)
  G_stat_ix2 = get_G_stat_ix(L, P2, aa_list, P_MSA2, kT)

  # Get the vectors of statistical energies
  # for each alignment position, summed over all 20 amino acids
  # A measure of how out of equilibrium each alignment position is
  G_stat_i1 = get_G_stat_i(L, P1, P_MSA1, kT)
  G_stat_i2 = get_G_stat_i(L, P2, P_MSA2, kT)
  result = perturb(aln, aa_list, P, P_MSA, L, N, A, kT)
  return result
  


def get_aa_counts(L, aln, aa_list):
  F = []  # A list of 20-long unit vectors of aa counts for each alignment column
  # For each column
  for i in range(L):
    col = aln.get_column(i)
    f = zeros(20)
    for a in col:
      if a in aa_list:  # If it is a valid amino acid
        k = aa_list.index(a)
        f[k] += 1
    F.append(f)
  return F



def get_P(L, aln, aa_list, F, N, A):
  P = []
  for j in range(L):
    p = zeros(20,float)
    for a in aa_list:
      x = aa_list.index(a)
      # Work in log space for these densities
      p[x] = r.dbinom(F[j][x], N, A[x], TRUE)
    P.append(p)
  return P


  
def get_H(L, aln, aa_list):
  H = zeros(20, float)
  for i in range(L):
    col = aln.get_column(i)
    for a in col:
      if a in aa_list:
        x = aa_list.index(a)  # Returns a unique index in aa_list
        H[x] += 1
  # Now modify H. Turn int from counts overall in the alignment, to expected number given
  # overall frequency in alignment
  # H /= L   # This is the same as H = (H * N) / (L * N) or freq * number of sequences
  # H[x] now has the expected number of observations of amino acid x in the hypothetical column
  # They need to be ints for binomial function. So use normal approximation to binomial
  for x in range(20):
    h = float(H[x]) / L
    H[x] = h  # The expected number of observations
  return H



def get_P_MSA(N, A, H):
  P_MSA = zeros(20, float)
  for x in range(20):
    #P_MSA[x] = r.dbinom(H[x], N, A[x], TRUE)
    # Approximate binomial with normal distribution
    # We know H[x] the expected number of counts for the amino acid
    # We also need mean and variance for the binomial, then we will turn this
    # to parameters for the normal (mean and SD)
    # Mean in binomial is np
    mu = N * A[x]
    sd = pow(N * A[x] * (1-A[x]),.5)
    P_MSA[x] = r.dnorm(H[x], mu, sd, TRUE)
  return P_MSA



def get_G_stat_ix(L, P, aa_list, P_MSA, kT):
  G_stat_ix = []
  for j in range(L):  # For each column, look at amino acid frequencies
    p = P[j]
    g = zeros(20,float)
    for x in range(20):
      a = aa_list[x]
      arg = p[x] - P_MSA[x]  # Since they are logs this difference is the log of the ratio
      g[x] = kT * arg
    G_stat_ix.append(g)
  return G_stat_ix



def get_G_stat_i(L, P, P_MSA, kT):
  G_stat_i = zeros(L, float)
  for i in range(L):
    sum = 0.0
    for x in range(20):
      sum += pow(P[i][x] - P_MSA[x], 2)
    G_stat_i[i] = kT * pow(sum, .5)
  return G_stat_i



def perturb(aln, aa_list, P, P_MSA, L, N, A, kT):
  I = range(aln.get_alignment_length())
  J = range(aln.get_alignment_length())

  result = zeros((len(I), len(J)), float)

  # For each column in the alignment
  for i in I:
    # For each other column in the alignment:
    for j in J:
      if i != j:
        sum = 0.0
	print
        print i, aln.get_column(i)
        print j, aln.get_column(j)
        # Generate all possible perturbations at site J
        # For each kind of amino acid
        for x in range(len(aa_list)):
          a = aa_list[x]
	  print "  %s" % a

          # Generate perturbation
          # Make another alignment, will be a subset with one amino acid held constant at some column
          subset = ClustalAlignment()   # Let's hope this doesn't result in unintended changes to aln  
          #subset._records = []  # Empty the list
          for s in range(len(aln._records)):
            seq = aln._records[s]
	    # Bug alert!
            if seq.seq[j] == a:  # If this sequence has amino acid x at position i
	      subset._records.append(seq)
          # Now get P and P_MSA for the perturbed alignment
          # Only if there is at least one sequence in the alignment though
          # Wouldn't make sense to calculate energy change for alignment with no sequences
	  if len(subset._records) > 0:
            F_delta = get_aa_counts(L, subset, aa_list)   # A list of 20-long unit vectors of aa counts for each alignment column
	    N_delta = len(subset._records)
            P_delta = get_P(L, subset, aa_list, F_delta, N_delta, A) 
            H_delta = get_H(L, subset, aa_list)
            P_MSA_delta = get_P_MSA(N_delta, A, H_delta)
	    a = (P_delta[i][x] - P_MSA_delta[x])
	    b = (P[i][x] - P_MSA[x])
            change = a - b 
            sum += change
            m = subset.get_column(i) 
            n = subset.get_column(j) 
	    print "    %s" % (m)
	    print "    %s" % (n)
	    print "    %f" % a
	    print "    %f" % b
	    print "    %f" % change
	print "    %f" % sum
        #result[i][j] = kT * pow(sum, .5)
        result[i][j] = sum
  return result



def get_aa_freq_from_db(file):
  aa_dict = {}
  num_aa = 0  # The total number of amino acid residues looked at
  for line in file.readlines():
    header_pat = search(r'^>.+$', line)
    data_pat = search(r'^(\w+)$', line)
    if data_pat:
      seq_frag = data_pat.group(1)
      for a in seq_frag:
        if a in aa_dict.keys():
	  aa_dict[a] += 1
	else:
	  aa_dict[a] = 1
	num_aa += 1
  # Divide all by number of amino acids in database
  for a in aa_dict.keys():
    aa_dict[a] /= float(num_aa)
  return aa_dict
  

  
################################################
# Some general matrix-number-crunching functions
################################################

def log_transform_matrix(A):
  X = zeros(A.shape, float)
  for i in range(len(A)):
    for j in range(len(A[i])):
      X[i][j] = -log(A[i][j])
  return X



def calc_matrix_correlation(A, B):  # Takes python Numeric matrix as argument
  # Turn each matrix into a unit vector
  X = zeros(product(A.shape), float)
  Y = zeros(product(B.shape), float)
  i = 0
  for row in A:
    for element in row:
      X[i] = element
      i += 1
  i = 0
  for row in B:
    for element in row:
      Y[i] = element
      i += 1
  r = pearsonr(X,Y)[0]
  return r


      
def matrix_to_unit_vector(A):
  X = zeros(product(A.shape), float)
  i = 0
  for row in A:
    for element in row:
      X[i] = element
      i += 1
  return X



def square_matrix_to_unit_vector(A):
  I = A.shape[0]
  J = A.shape[1]
  if I != J:
    print "ERROR: Your matrix is not square"
    sys.exit()
  X = zeros((I * (J-1) / 2), float) 
  x = 0
  for i in range(I):
    for j in range(J):
      if i > j:  # Only do the unique comparisons
        X[x] = A[i][j]
	x += 1
  return X



def calc_weighted_matrix_correlation(A, B, C):  
  # Takes 4 python Numeric matrices as arguments
  # All matrices are the same dimensions square
  # A and B are the ones we want to find the correlations between
  # w1 and w2 are the weights for A and B respectively
 
  # Check the matrix dimensions.
  for x in (A, B, C):
    for y in (A, B, C):
      if x.shape != y.shape:
        print "Matrices wrong shape"
	sys.exit()
  
  # Turn each matrix into a unit vector
  X = matrix2unit_vector(A)
  Y = matrix2unit_vector(B)
  w = matrix2unit_vector(C)

  N = len(X)  # And also is equal to the length of the other 3 vectors
      
  # Get the necessary means and standard deviations
  Mx = mean(X)
  My = mean(Y)
  Sx = sd(X)
  Sy = sd(Y)
  
  # Start the calculation of the correlation coefficient
  numerator = 0.0
  for i in range(N):
    arg = (w[i] * (X[i] - Mx) * (Y[i] - My)) / (Sx * Sy)
    numerator += arg
  r = numerator / N
  return r



def matrix2unit_vector(A):
  X = zeros(product(A.shape), float)
  i = 0
  for row in A:
    for element in row:
      X[i] = element
      i += 1
  return X 



def list2unit_vector(A):
  N = len(A)
  X = zeros(N, float)
  for i in range(N):
    X[i] = A[i]
  return X



def product(value_list):
  prod = 1
  for v in value_list:
    prod *= v
  return prod



def mean(X):
  N = float(len(X))
  sum = 0.0
  for x in X:
    sum += x
  mu = sum / N
  return mu



def sd(X):
  N = float(len(X))
  mu = mean(X)
  sum = 0.0
  for x in X:
    sum += pow(x - mu, 2)
  sigma = pow(sum / N, .5)  
  return sigma



def correlation(X,Y):
  if len(X) == len(Y):
    N = len(X)
  else:
    print "Unequal length vectors"
    sys.exit()
  Mx = mean(X)
  My = mean(Y)
  Sx = sd(X)
  Sy = sd(Y)
  numerator = 0.0
  for i in range(N):
    arg = ((X[i] - Mx) / Sx) * ((Y[i] - My) / Sy)
    numerator += arg
  r = numerator / N
  return r



def calc_z(r, N):
  from math import log
  z = (.5 * log((1 + r) / (1 - r))) / (pow( 1.0 / (N - 3), .5))
  return z


  
def calc_pval_from_z(z):
  z = fabs(z)  # Make it the absolute value
  pval = 1 - r.pnorm(z)
  return pval



def get_sander_weights(aln):  
  # Calculate the weight w used by Pazos et al., from Sander and Schneider Proteins 1991
  # Takes an alignment as a BioPython Alignment object
  # The length of the alignment
  L = aln.get_alignment_length()
  # The number of sequences
  N = len(aln._records)
  # An N x N matrix
  weights = zeros((N,N), float)
  # Get the summations of residue matches
  for i in range(L):
    sys.stdout.flush()
    col = aln.get_column(i)
    for k in range(N):
      for l in range(N):
        if k > l:
	  if col[k] == col[l]:
	    weights[l][k] += 1.0
	    weights[k][l] += 1.0
	elif k == l:
	  weights[k][l] += 1.0
  # Divide by L
  weights /= L
  weights = 1 - weights
  return weights



def get_empty_weights(aln):
  N = len(aln._records)
  # An N x N matrix
  weights = ones((N,N), float)
  return weights



def concatenate_alignments(aln1, aln2):
  joined_aln = ClustalAlignment()
  if len(aln1._records) != len(aln2._records):
    print "Different numbers of sequences in your two alignments"
    sys.exit()
  for k in range(len(aln1._records)):
    a = aln1._records[k].seq
    b = aln2._records[k].seq
    c = a + b
    d = aln1._records[k].description
    e = aln2._records[k].description
    f = d + '-' + e
    joined_aln._records.append(SeqRecord(c, description=f))
  return joined_aln


def combine_sander_weights(w1, w2):
  # Come up with a single weight for both alignments.
  # Position l,k is the weight for comparing sequence l to sequence k
  # This is fine for intra-protein correlation analysis. But what about
  # inter-protein? For now, just take the product of the weights.
  # Another way would be to concatenate the alignments and then get the
  # fraction mismatches over both.
  return w1 * w2



def running_average(X,Y,out):
  bins = make_bins(0.0, 80.0, 5.0)
  bin_data(bins, X, Y)
  A, B = get_running_average(bins)
  ra_cor = r.cor(A, B, method = "pearson")
  print "\tRunning average:\t%f" % ra_cor
  write_running_average(bins, out)



def make_bins(d_init, d_max, window):
  bins = {}
  i = d_init 
  while i < d_max:
    bins[(i, i + window)] = []  # List of correlation values
    i += 1.0
  return bins



def bin_data(bins, X, Y):
  if len(X) != len(Y):
    print "ERROR: vectors different lengths"
    sys.exit()
  N = len(X)
  for i in range(N):
    x = X[i]
    y = Y[i]
    for a, b in bins.keys():
      if x > a and x <= b:
        bins[(a,b)].append(y)
  return bins



def write_running_average(bins, file):
  sorted_keys = bins.keys()
  sorted_keys.sort()
  for a, b in sorted_keys:
    n = len(bins[(a,b)])
    if n > 0:
      m = mean(bins[(a,b)])
      file.write("%.1f\t%.1f\t%d\t%f\n" % (a, b, n, m))



def get_running_average(bins):
  min_points = 5
  A = []
  B = []
  for a, b in bins.keys():
    n = len(bins[(a,b)])
    if n >= min_points:
      m = mean(bins[(a,b)])
      A.append(a)
      B.append(m)
  return A, B



def fodor_test(scores,dist):
  dist_percentile = 50
  prob = dist_percentile / float(100)
  dist_threshold = scoreatpercentile(dist, dist_percentile)
  score_rank = 75
  
  data = []
  for i in range(len(scores)):
    data.append((scores[i], dist[i]))
  data.sort()     # Sorts on first element of tuple
  data.reverse()  # Low to high

  num_proximal = 0
  for i in range(score_rank):
    s = data[i][0]
    d = data[i][1]
    if d <= dist_threshold:
      num_proximal += 1
  pval = 1 - binom.cdf(num_proximal - 1, score_rank, prob) 
  print "%d\t%d\t%d\t" % (num_proximal, score_rank, dist_percentile),
  print "%f" % pval



def angstrom_test(scores,dist):
  dist_threshold = 15.0  # The number of angstroms
  score_ranks = (25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000)
  print
  print "Distance cutoff: %f" % dist_threshold
  print "Score ranks: ", score_ranks
  
  data = []
  for i in range(len(scores)):
    data.append((scores[i], dist[i]))
  data.sort()     # Sorts on first element of tuple
  data.reverse()  # Low to high
  N = len(data)
  for rank in score_ranks:
    if rank <= N:
      num_proximal = 0
      score_threshold = data[rank][0]
      k = 0
      m = 0
      n = 0
      for i in range(N):
        s = data[i][0]
        d = data[i][1]
        if s >= score_threshold and d <= dist_threshold:
          k += 1
        if s >= score_threshold:
          m += 1
        if d <= dist_threshold:
          n += 1
      pval = r.phyper(k-1, n, N-n, m, 0) 
      tp = k
      fn = n - k
      fp = m - k
      tn = N + k - m - n
      sens = tp / float(tp + fn)
      spec = tn / float(tn + fp)
      ppv = tp / float(tp + fp)
      f_m = m / float(N)
      f_n = n / float(N)
      exp = f_m * f_n * N
      obs = k
      enrichment = obs / exp
      print "%d\t%d\t%f\t%f\t%f" % (k, rank, dist_threshold, enrichment, pval)
      sys.stdout.flush()



def percentile_test_hypergeometric(scores,dist):
  significance_threshold = .05
  print
  dist_percentiles = [1, 2, 5, 10, 20, 50, 75, 99]
  score_percentiles = [99, 98, 95, 90, 80, 50, 25, 0]
  score_probs = array(score_percentiles, float)
  score_probs /= 100.0
  dist_probs = array(dist_percentiles, float)
  dist_probs /= 100.0
  score_thresholds = r.quantile(scores, score_probs, names=False) 
  print score_thresholds
  dist_thresholds = r.quantile(dist, dist_probs, names=False) 
  print dist_thresholds

  data = []
  for i in range(len(scores)):
    data.append((scores[i], dist[i]))
  N = len(data)

  for x in score_percentiles:
    for y in dist_percentiles:
      score_threshold = score_thresholds[score_percentiles.index(x)]
      dist_threshold = dist_thresholds[dist_percentiles.index(y)]
      k = 0
      m = 0
      n = 0
      for i in range(N):
        s = data[i][0]
        d = data[i][1]
        if s >= score_threshold and d <= dist_threshold:
          k += 1
	if s >= score_threshold:
	  m += 1
	if d <= dist_threshold:
	  n += 1
      pval = r.phyper(k-1, n, N-n, m, 0) 
      if pval <= significance_threshold:
        result = 1
      else:
        result = 0
      tp = k
      fn = n - k
      fp = m - k
      tn = N + k - m - n
      sens = tp / float(tp + fn)
      spec = tn / float(tn + fp)
      ppv = tp / float(tp + fp)
      f_m = m / float(N)
      f_n = n / float(N)
      exp = f_m * f_n * N
      obs = k
      enrichment = obs / exp
      print "%d\t%f\t%d\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%.30f" % (x, score_threshold, m, y, dist_threshold, n, k, N, sens, spec, ppv, enrichment, pval)
      sys.stdout.flush()



def calc_euclidean_distance(A, B):
  if len(A) == len(B):
    L = len(A)
  else:
    print "Unequal length vectors"
  sum = 0.0
  for i in range(L):
    sum += pow(A[i] - B[i], 2)
  return pow(sum, .5)



def minimize_vectors(A,B):  # Truncate the vectors to those that are nonzero for both positions
  if len(A) == len(B):
    L = len(A)
  else:
    print "Unequal length vectors"
  # Copy the vectors leaving out pairs of measurements that are zero for both
  # What is the number where both aren't zero? Use for making zeros vector
  M = []
  N = []
  for i in range(L):
    if A[i] != 0.0 or B[i] != 0.0:  # if either is non-zero
      # For coevolution matrices, this excludes the case when A[i] and B[i] are both zero.
      # That will only happen when we evaluate sequences that are identical to a pair in the
      # reference matrix
      M.append(A[i])
      N.append(B[i])
  X = array(M)
  Y = array(N)
  return X, Y



def distmat_correlation(protlist1, protlist2, ddict1, ddict2):
  A = []
  B = []
  for i in protlist1:
    for j in protlist1:
      A.append(ddict1[i][j])
  for i in protlist2:
    for j in protlist2:
      B.append(ddict2[i][j])
  r = correlation(A, B)
  return r
