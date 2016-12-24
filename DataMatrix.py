#! /usr/bin/env python
#-*- coding:utf-8 -*-

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class DataMatrix:
  def __init__(self, eps, filename):
    
    self.eps = eps
    self.data = pd.read_csv(filename, decimal=",", sep="\s+")

    self.headers = list(self.data.columns.values)

    for col in self.headers[1:]:
    
      #Turn the string values to float
      self.data[col] = pd.to_numeric(self.data[col], errors = 'coerce') 
    
      #Normalize data
      self.data[col] = (self.data[col] - self.data[col].mean())/(self.data[col].max()-self.data[col].min())

    self.nb_el = len(self.data)

    self.L = np.array([[0.0  for _ in range(self.nb_el)] for _ in range(self.nb_el)])

    for i in range(self.nb_el):
      for j in range(self.nb_el):
        self.L[i,j] = self.k(i,j)

  def k(self, i ,j):
    """Returns the gaussian similarity kernel between the elements at positions i and j"""

    x = self.data.iloc[i,1:]
    y = self.data.iloc[j,1:]

    return ((np.linalg.norm(x-y))**2)/self.eps


if __name__ == '__main__':
  
  test = DataMatrix(2, "data_small.csv")
  print test.k(2,3)

