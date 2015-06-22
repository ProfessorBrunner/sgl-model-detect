#! /usr/bin/env python
#@Author Sean McLaughlin

desc = '''
This is a helper function that is not called by the main modeling function. It converts the Ra, Dec of an object into pixel coordinates, and saves that file in the right format
In the directory. That file is the "centers" file used in the main modeler.
'''
import argparse
parser = argparse.ArgumentParser(description = desc)

parser.add_argument('dirname', metavar = 'dirname', type = str, help = 'Directory from which to read the .fits files and to write the output file.')
parser.add_argument('catalog', metavar = 'catalog', type = str, help = 'File where the informatoin regarding the .fits images is located. ')
parser.add_argument('ra_idx', metavar = 'ra_idx', type = int, default= 2, help = 'Index of the ra in the catalog along a line. Assumed to be followed by dec.', nargs = '?')
parser.add_argument('run_idx', metavar = 'run_idx', type = int, default = 14, help= 'Index of the run in the catalog along a line. Assumed to be followed by rerun, camcol, and field.', nargs = '?')

args = parser.parse_args()

raDecDict = {}
with open(args.catalog) as f:
    for line in f:
        line = line.strip()
        if line[0] == '#':
            continue
        splitLine = line.split(' ')
        #Files sometimes have undefined splits between entries. This removes them dilligently.

        try:
            while True:
                splitLine.remove(' ')
        except ValueError:
            pass

        ra, dec = splitLine[args.ra_idx:args.ra_idx+2]
        run, rerun, camcol, field = splitLine[args.run_idx:args.run_idx+4]
        for var in (ra, dec, run, rerun, camcol, field):
            var = var.strip()

        id = '-'.join(['0'*(6-len(run))+run, camcol, '0'*(4-len(field))+field])
        print id