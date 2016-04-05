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

md = np.random.rand(n)
ud = np.random.rand(n-1)
ld = np.random.rand(n-1)
rhs = np.random.rand(n)
w = np.random.rand(n,2)
v = np.random.rand(n,2)


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
sb.call(['./a.out',str(n)])

ans = np.fromfile('test_results.dat')

err= abs( (sol-ans[:n])/sol).max()
err1=abs( (sol1-ans[n:2*n])/sol1).max()
err2=abs( (sol2-ans[-n:])/sol2).max()

errt = np.array( [err,err1,err2])



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


else:
    print 'All tests sucessful'
    print '    Tri-diag: %.3e'%err
    print '    Tri-diag 1 weight: %.3e'%err1
    print '    Tri-diag 2 weights: %.3e'%err2





