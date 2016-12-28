#! /usr/bin/env python
#-*- coding:utf-8 -*-

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp

class DataMatrix:
  def __init__(self, eps, filename):
    
    self.eps = eps
    self.data = pd.read_csv(filename, decimal=",", sep="\s+")

    #Store column titles
    self.headers = list(self.data.columns.values)

    self.data = np.array(self.data.values)
    self.nb_el = len(self.data)

    
    for i, row  in enumerate(self.data[:,1:]) :
    
      #Turn the string values to float
      for j, elt in enumerate(row):  
        elt = float(elt)
      #Normalize data
        self.data[i,j+1] = (self.data[i,j+1] - self.data[:,j+1].mean())/(self.data[:,j+1].max()-self.data[:,j+1].min()) if self.data[:,j+1].max()-self.data[:,j+1].min() != 0 else 0 #(self.data[i,j] - self.data[:,j].mean())


    #Initialisation and making of L matrix
    self.L = np.array([[0.0  for _ in range(self.nb_el)] for _ in range(self.nb_el)])

    for i in range(self.nb_el):
      for j in range(self.nb_el):
        if i >= j :
          self.L[i,j] = self.k(i,j)
          self.L[j,i] = self.L[i,j]


    print self.L



  def k(self, i ,j):
    """Returns the gaussian similarity kernel between the elements at positions i and j"""
    x = self.data[i,1:]
    y = self.data[j,1:]
    return exp(((np.linalg.norm(x-y))**2)/self.eps)


if __name__ == '__main__':
  
  test = DataMatrix(2, "data_small.csv")
  

