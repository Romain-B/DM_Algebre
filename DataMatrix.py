#! /usr/bin/env python
#-*- coding:utf-8 -*-

import csv
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp

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
  def __init__(self, eps, filename):

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

    v,w = np.linalg.eig(self.Ms)
    le = len(v)
    v = np.sort(v)[::-1]

    #Initialisation de V et Lambda, Lambda contenant les valeurs propres de M 
    #et V les vecteurs propres.


    self.Lambda = np.array([[0.0 for _ in range(le)] for _ in range(le)])
    for i in range(le):
      self.Lambda[i,i]=v[i]
      sys.stdout.write(load_bar(i+1,le))

    self.V = w


    print "\nDone.\n"


    #Initialisation of M matrix :
    self.M = np.dot(D_mhalf,self.V).dot(self.Lambda).dot(self.V.transpose()).dot(D_half)



  def k(self, i ,j):
    """Returns the gaussian similarity kernel between the elements at positions i and j"""
    x = self.data[i,1:]
    y = self.data[j,1:]
    return exp(-((np.linalg.norm(x-y))**2)/self.eps)


if __name__ == '__main__':
  
  test = DataMatrix(0.5, "data_small.csv")
  print "\nM :", test.M
  # print "Valeurs propres de Ms :", test.Lambda
  # print "Vecteurs propres de Ms :", test.V
  #test = DataMatrix(2, "data.csv")

