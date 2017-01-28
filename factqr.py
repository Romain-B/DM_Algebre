#!/usr/bin/python2.6
# -*-coding:Latin-1 -*

from numpy import *

"""
module python pour factorisation QR avec méthodes de Gram-Schmidt et Householder
Algèbre linéaire et matricielle 3BIM INSA Lyon
2014 v20/11/2014
"""

def clgs(A):
    """factorisation QR Gram-Schmidt classique - instable
    
    Usage:
    Q, R = clgs(A)
    
    Argument:
    A: matrice mxn avec m>=n
    
    Retour:
    Q: matrice Q de la factorisation QR
    R: matrice R de la factorisation QR
    
    Note:
    Cette factorisation n'est pas stable numériquement
    Tire de: LN Trefethen & D Bau III, Numerical Linear Algebra, 1997 SIAM Philadelphia 
    """
    
    m,n = A.shape
    v = zeros((m,n))
    R = zeros((n,n))
    Q = zeros((m,n))

    for j in range(0,n):
        v[:,j] = A[:,j]
        for i in range(0,j):
            R[i,j] = dot( (Q[:,i].conjugate()).T, A[:,j] )
            v[:,j] = v[:,j] - R[i,j]*Q[:,i]
            
        R[j,j] = linalg.norm(v[:,j],2)
        Q[:,j] = v[:,j]/R[j,j]

    return (Q,R)

def mgs(A):
    """factorisation QR Gram-Schmidt modifiée - stable
    
    Usage:
    Q, R = mgs(A)
    
    Argument:
    A: matrice mxn avec m>=n
    
    Retour:
    Q: matrice Q de la factorisation QR
    R: matrice R de la factorisation QR
    
    Note:
    Cette factorisation est stable numériquement
    Tire de: LN Trefethen & D Bau III, Numerical Linear Algebra, 1997 SIAM Philadelphia 
    """
    
    m,n = A.shape
    v = zeros((m,n))
    R = zeros((n,n))
    Q = zeros((m,n))

    for i in range(0,n):
        v[:,i] = A[:,i]

    for i in range(0,n):
        R[i,i] = linalg.norm(v[:,i],2)
        Q[:,i] = v[:,i]/R[i,i]
        for j in range(i,n):
            R[i,j] = dot( (Q[:,i].conjugate()).T, v[:,j])
            v[:,j] = v[:,j] - R[i,j]*Q[:,i]

    return (Q,R)

def  householder(A): # computes a QR factorization upper triangular matrix R 
    """matrices W, R pour la factorisation QR'
    
    Usage: 
    W, R = householder(A)

    Argument:
    A: matrice mxn avec m>=n
    
    Retour:
    W: matrice de vecteurs de reflexions de householder
    R: matrice R de la factorisation QR

    Tire de: LN Trefethen & D Bau III, Numerical Linear Algebra, 1997 SIAM Philadelphia 
    """
    m,n = A.shape
    ident = eye(m)
    W = zeros((m,n))
    R = copy(A)

    for k in range(0,n):
        e1=ident[k:m,k]
        x = R[k:m,k]
        W[k:m,k]= ( int(x[0]>=0) - int(x[0]<0) )*linalg.norm(x,2)*e1+x
        W[k:m,k]=W[k:m,k]/linalg.norm(W[k:m,k],2)
        R[k:m][:,k:n]=R[k:m][:,k:n]-dot(outer(2*W[k:m,k],W[k:m,k]),R[k:m][:,k:n])

    return W, R

def formQ(W): 
    """Matrice Q de la factorisation QR
    
    Usage: 
    Q = formQ(W)
    
    Argument:
    W: matrice mxn de reflexions obtenue de la décomposition Householder
    
    Retour:
    Q: Matrice Q de la factorisation QR
    
    Note:
    La matrice W se calcule avec W, R = householder(A)
    """
    
    m,n = W.shape
    Q=zeros((m,m))

    # compute the product QI = Q, column by column
    for i in range(0,m):
        Q[i,i]=1
        for k in range(n-1,-1,-1): 
            Q[k:m,i]=Q[k:m,i]-2*dot(outer(W[k:m,k],W[k:m,k]),Q[k:m,i])

    return Q

def backsubs(R,b):
    """Solution d'un système linéaire - méthode substitution inverse
    
    Usage:
    x = backsubs(R,b)
    
    Arguments:
    R: matrice carrée mxm triangulaire supérieure inversible
    b: vecteur mx1
    
    Retour:
    x: solution de Rx = b
    """

    m, n = R.shape
    x = zeros((m,1))
    
    for j in range(m-1,-1,-1):
        xr = dot(R[j,(j+1):m],x[(j+1):m])
        x[j] = (b[j]-xr)/R[j,j]

    return x

def Qadjb(W,b): 
    """Produit matrice-vecteur Q adjoint fois b
    
    Usage:
    qb = Qadjb(W,b)
    
    Arguments:
    W: matrice mxn de reflexions obtenue de la décomposition Householder
    b: vecteur mx1
    
    Retour:
    Produit Q adjoint par b, où Q est la matrice de la factorisation QR
    
    Note:
    La matrice W se calcule avec W, R = householder(A)
    Tire de: LN Trefethen & D Bau III, Numerical Linear Algebra, 1997 SIAM Philadelphia 
    """

    m,n = W.shape
    qb = copy(b)

    for k in range(0,n):
        qb[k:m]=qb[k:m]-2*dot(outer(W[k:m,k],W[k:m,k]),qb[k:m])

    return qb



