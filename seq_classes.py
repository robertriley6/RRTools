####################################################
# VARIOUS FUNCTIONS FOR DEALING WITH BIOLOGICAL DATA
####################################################

import sys
import re
from re import findall, search, split
from math import log


###############
# THIS AND THAT
###############

TRUE = 1
FALSE = 0
ERROR = -1



#########
# CLASSES
#########

class Seq:
  def __init__(self, h, d):
    self.header = h
    self.data = d  # Usually a string
    self.id = ''   # Unless otherwise specified.  Usually put genome id here.


class Alignment:
  def __init__(self):
    self.sequences = []

  def write(self):
    incr = 60
    for s in self.sequences:
      print ">%s" % s.header
      a = 0
      L = len(s.data)
      while a < (L - incr):
        print "%s" % s.data[a:a+incr]
      	a += incr
      print "%s" % s.data[a:]         

  def write_string(self):
    str = ''
    for s in self.sequences:
      str += ">%s\n" % s.header
      str +=  "%s\n" % s.data
    return str
  
  def get_column(self, index):
    n = index - 1  # We usually think of 1 as the first number in a sequence, but this is not how arrays are indexed
    col = ''
    for s in self.sequences:
      col += s.data[n]
    return col 
  
  def length(self):
    L = len(self.sequences[0].data)
    for s in self.sequences:
      M = len(s.data)
      if M != L:
        print "You got a problem in your alignment"
        return ERROR
    return L

# A syntenic block
class SyntenicBlock:  # Could also extend the alignment class
  def __init__(self):
    self.ref_gen = ''
    self.sequences = []  # To be a list of Seq objects
    self.number = 0  # This references back to a number from a synteny block run.  For example in alignment_catalog.csv
    self.start = 0  # Start genome coord from H37Rv seq in alignment file
    self.end = 0    # End coord
  def ref_seq(self):
    for s in self.sequences:
      if search(self.ref_gen, s.header):
        return s
    

# And a related class extending Seq.  
# The reason for this is it has a backreference to the genome, contig, and interval 
# plus the header and sequence data
class BlockSeq(Seq):
  def __init__(self, h, d, genome, c, a, b):
    self.header = h
    self.data = d
    self.genome = g
    self.contig_num = c
    self.start = a  # w.r.t contig
    self.end = b    # same


# That .ptt file from NCBI
class PTT_annot:
  def __init__(self, start, end, strand, length, pid, gene, synonym, code, cog, product):
    self.start = start
    self.end = end
    self.strand = strand	
    self.length = length	
    self.pid = pid
    self.gene = gene	 
    self.synonym = synonym	
    self.code = code	
    self.cog = cog	
    self.product = product
  
class VistaAnnot:
  def __init__(self, strand, start, end, label):
    self.strand = strand
    self.start = start
    self.end = end
    self.label = label
    

class CatalogInterval:
  def __init__(self, aln_num, length, start, end):
    def process_id(aln_num):
      # Prepend some zeros to make a six digit number
      id = 'align'
      for i in (1, 2, 3, 4, 5, 6):
        if len(aln_num) == i:
          for j in range(6 - i):
            id += '0'
      id += aln_num
      return id
    self.id = process_id(aln_num)  # Convert '77' to 'align000077' etc. 
    self.length = length
    self.start = start
    self.end = end

class Cat:
  def __init__(self, id, start, end):
    self.id = id
    self.start = start
    self.end = end


class GFF_annot:
  def __init__(self, seqname, source, feature, start, end, score, strand, frame, attributes, comments):
    self.seqname = seqname
    self.source = source
    self.feature = feature
    self.start = start
    self.end = end
    self.score = score
    self.strand = strand	
    self.frame = frame
    self.attributes = attributes
    self.comments = comments
  def print_str(self):
    print "%s\t%s\t%s\t%d\t%d\t%s\t%1s\t%s\t%s" % \
          (self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attributes), 
    if self.comments != '':
      print "%s" % self.comments,
    print
  def write(self, f):
    f.write("%s\t%s\t%s\t%d\t%d\t%s\t%1s\t%s\t%s" % \
           (self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attributes))
    if self.comments != '':
      f.write(" %s" %  self.comments)
    f.write("\n")

# Class for the disjoint subsequences I_k in Tompa and Ruzzo paper
class Subseq:
  def __init__(self, left, right, score):
    self.left = left
    self.right = right
    self.score = score

