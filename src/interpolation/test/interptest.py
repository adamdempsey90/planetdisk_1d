#!/usr/bin/env python

import numpy as np
from scipy.interpolate import interp1d
from subprocess import call

nd = 1e4
nx = 1e3

tol = 1e-15

x = np.linspace(-15.,30.,nd)
y = np.zeros(x.shape)
indL = x<-1.0
indR = x>1.0
y[indL] = -(1./(x[indL]))**(4)
y[indR] = (1./(x[indR]))**(4)

r = (1 + x*.05);

x1 = np.exp(np.linspace(np.log(min(r)*1.02),np.log(.98*max(r)),nx))

func = interp1d(r,y)
y1 = func(x1)

with open("data.dat","w") as f:
    for xd,yd in zip(r,y):
        f.write('%.16f\t%.16f\n'%(xd,yd))

with open("xvals.dat","w") as f:
    for xd in x1:
        f.write('%.16f\n'%xd)

call(['./compile'])
call(['./interptest','%d'%nd,'%d'%nx])

yt = np.loadtxt('yvals_lint.dat')

err = np.sqrt(sum((y1-yt)*(y1-yt)))
print 'Total L2 error = %.4e' % err
if err <= tol:
    print 'Test PASSED for linear interpolation'
else:
    print 'Test FAILED for linear interpolation'


yt = np.loadtxt('yvals_splint.dat')

func = interp1d(r,y,kind='cubic')
y1 = func(x1)
err = np.sqrt(sum((y1-yt)*(y1-yt)))
print 'Total L2 error = %.4e' % err
if err <= tol:
    print 'Test PASSED for cubic spline interpolation'
else:
    print 'Test FAILED for cubic spline interpolation'
