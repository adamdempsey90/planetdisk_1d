#!/usr/bin/env python

import numpy as np
import subprocess as sb

import sys

if len(sys.argv) != 2:
    n = 10
else:
    n = int(sys.argv[1])

tol = 1e-10
print 'Testing with %d points'%n
print 'Using a tolerance of %e'%tol
print 'Copying over trisolve.c from src'
sb.call(['cp','../trisolve.c','.'])
sb.call(['cp','../matrix_ops.c','.'])

md = np.random.uniform(-2,2,size=n)
ud = np.random.uniform(-2,2,size=n-1)
ld = np.random.uniform(-2,2,size=n-1)
rhs = np.random.uniform(-2,2,size=n)
fm = np.random.uniform(-2,2,size=n)
wz = np.zeros((n,2))
vz = np.zeros((n,2))
w = np.random.uniform(-2,2,size=(n,2))
v = np.random.uniform(-2,2,size=(n,2))
a = np.random.uniform(-2,2)
b = np.random.uniform(-2,2)

mat = np.diag(md,0) + np.diag(ld,-1) + np.diag(ud,1)
sol = np.linalg.solve(mat,rhs)
sol3 = b*np.dot(mat,rhs) + a*fm
mat += np.outer(v[:,0],w[:,0])
sol1 = np.linalg.solve(mat,rhs)
mat += np.outer(v[:,1],w[:,1])
sol2 = np.linalg.solve(mat,rhs)
sol4 = b*np.dot(mat,rhs) + a*fm


with open('test_data.dat','wb') as f:
    md.tofile(f)
    ld.tofile(f)
    ud.tofile(f)
    rhs.tofile(f)
    fm.tofile(f)
    w.transpose().tofile(f)
    v.transpose().tofile(f)
    np.array([a]).tofile(f)
    np.array([b]).tofile(f)

sb.call(['./compile'])
sb.call(['./tritest',str(n)])
ans = np.fromfile('test_results.dat')


err= abs( (sol-ans[:n])/sol).max()
err1=abs( (sol1-ans[n:2*n])/sol1).max()
err2=abs( (sol2-ans[2*n:3*n])/sol2).max()
err3=abs( (sol3-ans[3*n:4*n])/sol3).max()
err4=abs( (sol4-ans[4*n:5*n])/sol4).max()


errt = np.array( [err,err1,err2,err3,err4])



if np.any(errt > tol):
    print 'Some tests failed'

    if err > tol:
        print '    Tri-diag no weights failed: %.3e'%err
        if n < 100:
            for s,e in zip(sol,ans[:n]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag no weights passed'
    if err1 > tol:
        print '    Tri-diag 1 weight failed: %.3e'%err1
        if n < 100:
            for s,e in zip(sol1,ans[n:2*n]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag 1 weight passed'
    if err2 > tol:
        print '    Tri-diag 2 weights failed: %.3e'%err2
        if n < 100:
            for s,e in zip(sol2,ans[2*n:3*n]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag 2 weights passed'
    if err3 > tol:
        print '    Matvc 0  weights failed: %.3e'%err3
        if n < 100:
            for s,e in zip(sol3,ans[3*n:4*n]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag 2 weights passed'
    if err4 > tol:
        print '    Matvec 2  weights failed: %.3e'%err4
        if n < 100:
            for s,e in zip(sol3,ans[4*n:5*n]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag 2 weights passed'


else:
    print 'All tests sucessful'
    print '    Tri-diag: %.3e'%err
    print '    Tri-diag 1 weight: %.3e'%err1
    print '    Tri-diag 2 weights: %.3e'%err2
    print '    Matvec 0 weights: %.3e'%err3
    print '    Matvec 2 weights: %.3e'%err4





