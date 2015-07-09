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
parser.add_argument('centerFilename', metavar = 'centerFilename', type = str, help ='File to write the output to.')
parser.add_argument('ra_idx', metavar = 'ra_idx', type = int, default= 2, help = 'Index of the ra in the catalog along a line. Assumed to be followed by dec. Default is 2.', nargs = '?')
parser.add_argument('run_idx', metavar = 'run_idx', type = int, default = 14, help= 'Index of the run in the catalog along a line. Assumed to be followed by rerun, camcol, and field. Default is 14.', nargs = '?')

args = parser.parse_args()

raDecDict = {}

from os import listdir
from fnmatch import fnmatch
import astropy.io.fits as fits
import astropy.wcs as wcs

filenames = listdir(args.dirname)
with open(args.catalog) as f:
    for line in f:
        line = line.strip()
        if line[0] == '#':
            continue
        splitLine = line.split()

        #Grab important atriubtes
        ra, dec = splitLine[args.ra_idx:args.ra_idx+2]
        run, rerun, camcol, field = splitLine[args.run_idx:args.run_idx+4]
        for var in (ra, dec, run, rerun, camcol, field):
            var = var.strip()

        id = '-'.join(['0'*(6-len(run))+run, camcol, '0'*(4-len(field))+field])
        if id in raDecDict: #alrady got this same image
            continue

        #identify the file that matches it
        for file in filenames:
            if fnmatch(file, '*'+id+'.fits'):
                raDecDict[file] = (id, ra,dec)
                break

toWrite = []
for filename, (id, ra, dec) in raDecDict.iteritems():
    hdulist = fits.open(args.dirname+filename)
    w = wcs.WCS(hdulist[0].header, hdulist)
    hdulist.close()
    #Convert Ra, Dec to pixel
    x,y = w.wcs_world2pix(float(ra), float(dec), 1)
    toWrite.append(' '.join([id, '%.3f'%x, '%.3f'%y]))

with open(args.centerFilename, 'w') as f:
    f.write('\n'.join(toWrite))