#! /usr/bin/env python
#-*- coding:utf-8 -*-

import csv
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp, sqrt
from mpl_toolkits.mplot3d import Axes3D


def load_bar(now, total):
  out = "["
  prog = int(((now*1.0)/(total*1.0)*100))

  for _ in range(prog):
    out += "#"
  for _ in range(100-prog):
    out+= ' '
  out += "] "+str(prog)+"%" 
  return "\r"+out



class DataMatrix:
  def __init__(self, eps, filename, k=2, t=5):

    print "\nLoading Dataset..."
    
    self.eps = eps
    self.data = pd.read_csv(filename, decimal=",", sep="\s+")
    print "Done.\n"

    #Store column titles
    self.headers = list(self.data.columns.values)

    self.data = np.array(self.data.values)
    self.nb_el = len(self.data)

    print "\nNormalizing data..."
    for i, row  in enumerate(self.data[:,1:]) :
      #Turn the string values to float
      for j, elt in enumerate(row):  
        elt = float(elt)
      #Normalize data
        self.data[i,j+1] = (self.data[i,j+1] - self.data[:,j+1].mean())/(self.data[:,j+1].max()-self.data[:,j+1].min()) if self.data[:,j+1].max()-self.data[:,j+1].min() != 0 else 0 #(self.data[i,j] - self.data[:,j].mean())
      sys.stdout.write(load_bar(i+1, self.nb_el))
      
    print "\nDone.\n"

    #Initialisation and making of L matrix
    self.L = np.array([[0.0  for _ in range(self.nb_el)] for _ in range(self.nb_el)])

    print "\nLoading L matrix...:\n"
    for i in range(self.nb_el):
      for j in range(self.nb_el):
        if i >= j :
          self.L[i,j] = self.k(i,j)
          self.L[j,i] = self.L[i,j]
      sys.stdout.write(load_bar(i+1,self.nb_el))
    print "\nDone.\n"

    #Initialisation and making of D matrix
    self.D = np.array([[0.0  for _ in range(self.nb_el)] for _ in range(self.nb_el)])


    print "Loading D matrix...\n"
    for i in range(self.nb_el):
      self.D[i,i] = self.L[i].sum()
      sys.stdout.write(load_bar(i+1,self.nb_el))
    print "\nDone.\n"

    #Iitialisation and making of Ms matrix with Ms = D^-1/2 L D^-1/2

    print "Loading Ms matrix...\n"

    D_half = np.matrix(np.zeros(self.D.shape))
    D_mhalf = np.matrix(np.zeros(self.D.shape))

    np.fill_diagonal(D_half, 1/ (self.D.diagonal()**0.5)) 
    np.fill_diagonal(D_mhalf, 1/ (self.D.diagonal()**-0.5))

    self.Ms = np.dot(D_mhalf,self.L).dot(D_mhalf)

      #Initialisation with M = D^-1 L
      #self.M = np.linalg.inv(self.D).dot(self.L)
    print "\nDone.\n"


    print "Computing Eigenvalues of Ms...\n"

    w,v = np.linalg.eig(self.Ms)
    le = len(v)
    #v = np.swapaxes(v,0,1)

    print "w,v\n",w,"\n",v

    #Sort in >>
    
    sor = w.argsort()[::-1]
    v = v[sor]

    w = w[sor]#np.sort(w)[::-1]

    #Initialisation de V et Lambda, Lambda contenant les valeurs propres de M 
    #et V les vecteurs propres.


    self.Lambda = np.array([[0.0 for _ in range(le)] for _ in range(le)])
    for i in range(le):
      self.Lambda[i,i]=w[i]
      sys.stdout.write(load_bar(i+1,le))

    self.V = v


    print "\nDone.\n"


    #Initialisation of M matrix :
    self.M = np.dot(D_mhalf,self.V).dot(self.Lambda).dot(self.V.transpose()).dot(D_half)

    self.Ksi = np.array([[0.0 for _ in range(self.nb_el)] for _ in range(self.nb_el)])
    self.Phi = np.array([[0.0 for _ in range(self.nb_el)] for _ in range(self.nb_el)])

    
    print "Loading Phi and Ksi matrixes...\n"
    #Save the v transpose to save computing time

    Vt = self.V.transpose()

    for j in range(self.nb_el):
      self.Ksi[j] = self.V[j]/sqrt(self.D[j,j])
      self.Phi[j] = Vt[j]*sqrt(self.D[j,j])
      sys.stdout.write(load_bar(j+1,self.nb_el))

    print "\nDone.\n"

    self.phi_chelou = []
    

    for i in range(k):
      self.phi_chelou.append(w[i+1]**t * self.Ksi[i+1])


  def k(self, i ,j):
    """Returns the gaussian similarity kernel between the elements at positions i and j"""
    x = self.data[i,1:]
    y = self.data[j,1:]
    return exp(-((np.linalg.norm(x-y))**2)/self.eps)


def randrange(n, vmin, vmax):
    '''
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    '''
    return (vmax - vmin)*np.random.rand(n) + vmin





if __name__ == '__main__':
  
  test = DataMatrix(0.1, "data_small.csv", k=3, t=20)
  print "\nD : \n", test.D
  print "\nMs :", test.Ms
  #print "\nPhi.Ksi (should be 1) : ", np.dot(test.Phi, test.Ksi)
  print "Valeurs propres de Ms :\n", test.Lambda
  print test.phi_chelou
  # print "Vecteurs propres de Ms :", test.V
  #test = DataMatrix(2, "data.csv")



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  n = 100

  # For each set of style and range settings, plot n random points in the box
  # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
  for c, m in [('r', 'o'), ('b', '^')]:
    xs = test.phi_chelou[0]
    ys = test.phi_chelou[1]
    zs = test.phi_chelou[2]
    ax.scatter(xs, ys, zs, c=c, marker=m)

  ax.set_xlabel('X Label')
  ax.set_ylabel('Y Label')
  ax.set_zlabel('Z Label')

  plt.show()