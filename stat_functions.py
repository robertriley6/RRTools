import sys
from scipy import exp
from scipy.special import gamma
from re import search, split



########################################
# MATHEMATICAL AND STATISTICAL FUNCTIONS
########################################



def fishers_binary(A,B):
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
  p = choose(m,k) * choose(N-m, n-k) / float(choose(N,n))
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
    #print "  %f" % p
    print x, p
    sum += p
    #print "%d\t%f" % (x, p)
    #print " %f" % sum
  return sum

def choose(n, k):
  a = fact(n)
  b = fact(n - k)
  c = fact(k)
  return a / float(b * c)

def fact(n):
  # Calculate n!
  if n > 170:  # The max on this (or most) machines
    a = pow(2 * 3.14159 * n, .5) * pow(n / exp(1), n)
  else:
    #a =  r.gamma(n+1)
    a =  gamma(n+1)
  #print n, a
  return a

