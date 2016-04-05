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

md = np.random.uniform(-2,2,size=n)
ud = np.random.uniform(-2,2,size=n-1)
ld = np.random.uniform(-2,2,size=n-1)
rhs = np.random.uniform(-2,2,size=n)
wz = np.zeros((n,2))
vz = np.zeros((n,2))
w = np.random.uniform(-2,2,size=(n,2))
v = np.random.uniform(-2,2,size=(n,2))


mat = np.diag(md,0) + np.diag(ld,-1) + np.diag(ud,1)
sol = np.linalg.solve(mat,rhs)
mat += np.outer(v[:,0],w[:,0])
sol1 = np.linalg.solve(mat,rhs)
mat += np.outer(v[:,1],w[:,1])
sol2 = np.linalg.solve(mat,rhs)

with open('test_data.dat','wb') as f:
    md.tofile(f)
    ld.tofile(f)
    ud.tofile(f)
    rhs.tofile(f)
    w.transpose().tofile(f)
    v.transpose().tofile(f)

sb.call(['./compile'])
sb.call(['./tritest',str(n)])
ans = np.fromfile('test_results.dat')

with open('test_data.dat','wb') as f:
    md.tofile(f)
    ld.tofile(f)
    ud.tofile(f)
    rhs.tofile(f)
    wz.transpose().tofile(f)
    vz.transpose().tofile(f)

sb.call(['./tritest',str(n)])

mat = np.diag(md,0) + np.diag(ld,-1) + np.diag(ud,1)
solz = np.linalg.solve(mat,rhs)
mat += np.outer(vz[:,0],wz[:,0])
solz1 = np.linalg.solve(mat,rhs)
mat += np.outer(vz[:,1],wz[:,1])
solz2 = np.linalg.solve(mat,rhs)
ansz = np.fromfile('test_results.dat')

err= abs( (sol-ans[:n])/sol).max()
err1=abs( (sol1-ans[n:2*n])/sol1).max()
err2=abs( (sol2-ans[-n:])/sol2).max()

errz= abs( (solz-ansz[:n])/solz).max()
errz1=abs( (solz1-ansz[n:2*n])/solz1).max()
errz2=abs( (solz2-ansz[-n:])/solz2).max()

errt = np.array( [err,err1,err2,errz,errz1,errz2])



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
            for s,e in zip(sol2,ans[-n:]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag 2 weights passed'


    if errz > tol:
        print '    Tri-diag no zero weights failed: %.3e'%errz
        if n < 100:
            for s,e in zip(solz,ansz[:n]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag no zero weights passed'
    if errz1 > tol:
        print '    Tri-diag 1 zero weight failed: %.3e'%errz1
        if n < 100:
            for s,e in zip(solz1,ansz[n:2*n]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag 1 zero weight passed'
    if errz2 > tol:
        print '    Tri-diag 2 zero weights failed: %.3e'%errz2
        if n < 100:
            for s,e in zip(solz2,ansz[-n:]):
                print '    %lg    %lg'%(s,e)
    else:
        print '    Tri-diag 2 zero weights passed'
else:
    print 'All tests sucessful'
    print '    Tri-diag: %.3e'%err
    print '    Tri-diag 1 weight: %.3e'%err1
    print '    Tri-diag 2 weights: %.3e'%err2
    print '    Tri-diag zero weights: %.3e'%errz
    print '    Tri-diag 1 zero weight: %.3e'%errz1
    print '    Tri-diag 2 zero weights: %.3e'%errz2





