#! /usr/bin/env python
#@Author Sean McLaughlin
desc = '''
This program models the light profile of a galaxy in a fits image, subtracts out the model and inspects the residuals for evidence of lensing.

This is the "main" module, that imports and runs the other modules in this directory
'''
from cropImage import cropImage
from findCenter import findCenter
from mcmcFit import mcmcFit
from residualID import residualID
import argparse
parser = argparse.ArgumentParser(description = desc)
#the number of Gaussians to use? (I think I may want to use a chi2 test to find the best one.)
#an option to change the number of walkers, to increase precision?
#I need a series of optional arguments to have save checkpoints along the process.
#Among those I need:
#Store the plot of the ID'd residuals
parser.add_argument('filename', metavar = 'filename', type = str, help = 'Either a fits filename or a directory of files.')
parser.add_argument('outputdir', metavar = 'output', type = str, help = 'Location to store the programs output')
parser.add_argument('centers', metavar = 'centers', type = str, help = 'Either a filename or a comma separate pair of coordinates for x,y. Default is to use findCenter.', default = None, nargs = '?')
parser.add_argument('--cutout', dest = 'cutout', action = 'store_const', const = True, default = False,\
                     help = 'Save a .png of the original cutout to file.')
parser.add_argument('--chain', dest = 'chain', action = 'store_const', const = True, default = False,\
                     help = 'Store the markov chain in the output directory. OCCUPIES A LOT OF FILESPACE.')
parser.add_argument('--triangle', dest = 'triangle', action = 'store_const', const = True, default = False,\
                     help = 'Save a .png of a triangle plot of the chain.')
parser.add_argument('--subtraction', dest = 'subtraction', action = 'store_const', const = True, default = False,\
                     help = 'Store a .png of the model subtraction from the original image.')
parser.add_argument('--residuals', dest = 'residuals', action = 'store_const', const = True, default = False,\
                     help = 'Save a .png of the detected residuals')

args = parser.parse_args()

