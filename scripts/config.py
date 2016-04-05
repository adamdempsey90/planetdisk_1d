#!/usr/bin/env/python

from sys import argv
from subprocess import call


def create_defines_file(optfile='../setups/params.opt',sourcedir='../src/'):
    if sourcedir[-1] != '/':
        sourcedir += '/'
    with open(optfile,'r') as f:
        temp = [x.split('+') for x in f.readlines()]
        defs=[]
        for x in temp:
            if x[0] == '' and '#' not in x:
                def_str = x[-1].split('\n')[0]
                defs.append(def_str)
    if 'NONLOCAL' in defs:
        call(['cp',sourcedir+'integration/crank_nicholson_step_nl.c',sourcedir+'integration.c'])
        if 'SHOCK' in defs:
            call(['cp',sourcedir+'torques/dTr_shocks.c',sourcedir+'dTr.c'])
        else:
            call(['cp',sourcedir+'torques/dTr_dep.c',sourcedir+'dTr.c'])
    else:
        call(['cp',sourcedir+'integration/crank_nicholson_step.c',sourcedir+'integration.c'])
        if 'GAUSSIAN' in defs:
            call(['cp',sourcedir+'torques/dTr_gaussian.c',sourcedir+'dTr.c'])
        else:
            call(['cp',sourcedir+'torques/dTr_linear.c',sourcedir+'dTr.c'])

    with open(sourcedir + 'defines.h','w') as g:
		if defs != []:
			for x in defs:
				g.write('#define ' + x + '\n')
		else:
			g.write('\n')
    print 'Created the defines.h file:'
    call(['cat',sourcedir+'defines.h'])
    return defs



if __name__ == "__main__":

	if len(argv) == 2:
		if str(argv[1])[-4:] == '.opt':
			optfile = str(argv[1])
			sourcedir = 'src/'
		else:
			sourcedir = str(argv[1])
			if sourcedir[-1] != '/':
				sourcedir += '/'
			optfile = 'params.opt'

	elif len(argv) == 3:

		if str(argv[1])[-4:] == '.opt':
			optfile = str(argv[1])
			sourcedir = str(argv[2])
		else:
			sourcedir = str(argv[1])
			if str(argv[2])[-4:] == '.opt':
				optfile = str(argv[2])
			else:
				print 'Could not find options file. Using default.'
				optfile = 'params.opt'

	else:
		print 'Using default options file and source directory'
		optfile = 'setups/params.opt'
		sourcedir = 'src/'

	if sourcedir[-1] != '/':
		sourcedir += '/'

	print 'Using Source Directory as: %s\nUsing Options File: %s' %(sourcedir,optfile)

	defs = create_defines_file(optfile,sourcedir)
