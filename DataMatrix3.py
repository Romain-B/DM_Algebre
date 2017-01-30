#! /usr/bin/env python
#-*- coding:utf-8 -*-

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import exp, sqrt


class DataMatrix:
  def __init__(self, eps, filename, k=2, t=5):
    
    self.eps = eps
    
#1. Import des données et normalisation   
#--------------
    self.data = pd.read_csv(filename, decimal=",", sep="\s+")

    #Stockage des en-têtes
    self.headers = list(self.data.columns.values)
    
    #Stockage des valeurs
    self.data = np.array(self.data.values)
    self.nb_el = len(self.data)

    
    for i, row  in enumerate(self.data[:,1:]) :
      for j, elt in enumerate(row):
        
        #Transformation des objets en float
        elt = float(elt)
        
        #Normalisation des données
        self.data[i,j+1] = (self.data[i,j+1] - self.data[:,j+1].mean())/(self.data[:,j+1].max()-self.data[:,j+1].min()) if self.data[:,j+1].max()-self.data[:,j+1].min() != 0 else 0 #(self.data[i,j] - self.data[:,j].mean())


#2. Construction de la matrice L
#--------------
    self.L = np.array([[0.0  for _ in range(self.nb_el)] for _ in range(self.nb_el)])

    for i in range(self.nb_el):
      for j in range(self.nb_el):
        if i >= j :
          self.L[i,j] = self.k(i,j)
          self.L[j,i] = self.L[i,j]

#3. Construction de la matrice D
#--------------
    self.D = np.array([[0.0  for _ in range(self.nb_el)] for _ in range(self.nb_el)])

    for i in range(self.nb_el):
      self.D[i,i] = self.L[i].sum()


#4. Construction de la matrice M avec M = D^-1/2 L D^1/2
#--------------
    D_half = np.matrix(np.zeros(self.D.shape))
    D_mhalf = np.matrix(np.zeros(self.D.shape))

    np.fill_diagonal(D_half,  (self.D.diagonal()**0.5))    # -> D^1/2
    np.fill_diagonal(D_mhalf, (self.D.diagonal()**-0.5))  # -> D^-1/2


    self.Ms = np.dot(D_mhalf,self.L).dot(D_mhalf)
    
    w,v = np.linalg.eig(self.Ms)
    le = len(v)
    
    sor = w.argsort()[::-1]
    v = v[sor]

    w = w[sor]
    
    self.Lambda = np.array([[0.0 for _ in range(le)] for _ in range(le)])
    for i in range(le):
      self.Lambda[i,i]=w[i]

    self.V = v

      #Construction de M si M = D^-1 L
      #self.M = np.linalg.inv(self.D).dot(self.L)
    
    #self.M = np.dot(D_mhalf,self.V).dot(self.Lambda).dot(self.V.transpose()).dot(D_half)

    self.Ksi = np.array([[0.0 for _ in range(self.nb_el)] for _ in range(self.nb_el)])
    self.Phi = np.array([[0.0 for _ in range(self.nb_el)] for _ in range(self.nb_el)])

    
    print "Loading Phi and Ksi matrixes...\n"
    #Save the v transpose to save computing time

    Vt = self.V.transpose()

    for j in range(self.nb_el):
      self.Ksi[j] = self.V[j]/sqrt(self.D[j,j])
      self.Phi[j] = Vt[j]*sqrt(self.D[j,j])
        
    self.phi_chelou = []
    

    for i in range(k):
      self.phi_chelou.append(w[i+1]**t * self.Ksi[i+1])

    

#5. Définition de la fonction noyau de similarité
#--------------
    
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


  datatest = DataMatrix(1,"data.csv", k=3, t=10)
  print "\nL :", datatest.L
  print "\nD :", datatest.D

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  n = 100

  # For each set of style and range settings, plot n random points in the box
  # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
  for c, m in [('r', 'o'), ('b', '^')]:
    xs = datatest.phi_chelou[2]
    ys = datatest.phi_chelou[1]
    zs = datatest.phi_chelou[0]
    ax.scatter(xs, ys, zs, c=c, marker=m)

  ax.set_xlabel('X Label')
  ax.set_ylabel('Y Label')
  ax.set_zlabel('Z Label')

  plt.show()

