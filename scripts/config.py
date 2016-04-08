#!/usr/bin/env/python

from sys import argv
from subprocess import call


def create_defines_file(optfile='../setups/params.opt',sourcedir='../src/',extras=None):
    if sourcedir[-1] != '/':
        sourcedir += '/'
    with open(optfile,'r') as f:
        temp = [x.split('+') for x in f.readlines()]
        defs=[]
        for x in temp:
            if x[0] == '' and '#' not in x:
                def_str = x[-1].split('\n')[0]
                defs.append(def_str)
#    if 'NONLOCAL' in defs:
#        call(['cp',sourcedir+'integration/crank_nicholson_step_nl.c',sourcedir+'integration.c'])
#        if 'SHOCK' in defs:
#            call(['cp',sourcedir+'torques/dTr_shocks.c',sourcedir+'dTr.c'])
#        else:
#            call(['cp',sourcedir+'torques/dTr_dep.c',sourcedir+'dTr.c'])
#    else:
#        call(['cp',sourcedir+'integration/crank_nicholson_step.c',sourcedir+'integration.c'])
#        if 'GAUSSIAN' in defs:
#            call(['cp',sourcedir+'torques/dTr_gaussian.c',sourcedir+'dTr.c'])
#        else:
#            call(['cp',sourcedir+'torques/dTr_linear.c',sourcedir+'dTr.c'])
    if extras is not None:
        defs += extras

    if 'OPENMP' in defs:
        add_library('OPENMP')
    else:
        rm_library('OPENMP')
    with open(sourcedir + 'defines.h','w') as g:
		if defs != []:
			for x in defs:
				g.write('#define ' + x + '\n')
		else:
			g.write('\n')
    print 'Created the defines.h file:'
    call(['cat',sourcedir+'defines.h'])
    return defs

def add_library(lib):
    print 'Adding %s to Makefile'%lib
    if lib.lower() == 'openmp':
        with open('Makefile','r') as f:
            lines = f.readlines()

        for i,line in enumerate(lines):
            if 'LDFLAGS=' in line:
                if 'OPENMP' not in line:
                    lines[i] = line.strip() + ' $(OPENMPLIB)\n'
            if 'CFLAGS=' in line:
                if 'OPENMP' not in line:
                    lines[i] = line.strip() + ' $(OPENMPFLAG)\n'
        with open('Makefile','w') as f:
            f.write(''.join(lines))
def rm_library(lib):
    print 'Removing %s to Makefile'%lib
    if lib.lower() == 'openmp':
        with open('Makefile','r') as f:
            lines = f.readlines()

        for i,line in enumerate(lines):
            if 'LDFLAGS' in line:
                if 'OPENMP' in line:
                    lines[i] = ' '.join(line.split('$(OPENMPLIB)'))
            if 'CFLAGS' in line:
                if 'OPENMP' in line:
                    lines[i] = ' '.join(line.split('$(OPENMPFLAG)'))
        with open('Makefile','w') as f:
            f.write(''.join(lines))


if __name__ == "__main__":

    sourcedir = 'src/'
    optfile = 'setups/params.opt'
    parfile = 'setups/params.par'
    extras = None
    for arg in argv[1:]:
        if '.opt' in arg:
            optfile = argv
        if '.par' in arg:
            parfile = arg
        if 'src' in arg:
            sourcedir = arg
            if  sourcedir[-1] != '/':
                sourcedir += '/'
        else:
            if extras is None:
                extras = [arg.upper()]
            else:
                extras.append(arg)
    print 'Using Source Directory as: %s\nUsing Options File: %s' %(sourcedir,optfile)
    if extras is not None:
        print 'Found these command line defines', extras
    defs = create_defines_file(optfile,sourcedir,extras)
