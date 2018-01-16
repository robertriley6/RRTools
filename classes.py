import sys

class pfam_hit:
  def __init__(self, prot_id, org, pf_id, pf_acc, pf_order, pf_count, start, end, score, eval):
    self.prot_id = prot_id
    self.org = org
    self.pf_id = pf_id
    self.pf_acc = pf_acc
    self.pf_order = pf_order
    self.pf_count = pf_count
    self.start = start
    self.end = end
    self.score = score
    self.eval = eval

