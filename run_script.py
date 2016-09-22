#!/usr/bin/env/python
from subprocess import call

xd_vals = [1.0,2.0,3.0,4.0,5.0]

fnames = ['results_xd%d.hdf5'%x for x in xd_vals]

with open('setups/params.in','r') as f:
    lines = f.readlines()

for x,fn in zip(xd_vals,fnames):
    print x,fn
    for i,line in enumerate(lines):
        if 'xd' in line:
            line = line.split('=')
            line[-1] = ' %.1f\n'%x
            lines[i] = '='.join(line)
            print lines[i]
        if 'outputname' in line:
            line = line.split('=')
            line[-1] = fn + '\n'
            lines[i] = '='.join(line)
            print lines[i]

    with open('setups/params.in.xd%d'%x,'w') as f:
        f.write(''.join(lines))

