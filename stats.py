import sys
from math import *

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

def corr(X,Y):
  if len(X) == len(Y):
    N = len(X)
  else:
    print "Unequal length vectors"
    print X
    print Y
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

def z_score(x, X):
  z = (x - mean(X)) / sd(X) 
  return z

def hamming(X, Y):
  L = len(X)
  if len(Y) != L:
    print "ERROR: vector lengths different: %d, %d" % (L, len(Y))
  else:
    h = 0
    for i in range(L):
      if X[i] != Y[i]:
        h += 1
  return h


def median(X):
  Y = sorted(X)

  if len(Y) % 2 == 1:
    return Y[(len(Y)+1)/2-1]
  else:
    lower = Y[len(Y)/2-1]
    upper = Y[len(Y)/2]

    return (float(lower + upper)) / 2  

def euclidean(P, Q):
  return sqrt(sum([pow(Q[i] - P[i], 2) for i in range(len(P))]))


def centroid(data):
  L = len(data[0])
  M = len(data)
  C = [0.0 for i in range(L)]
  for i in range(L):
    X = [x[i] for x in data]
    C[i] = mean(X)
  return C


