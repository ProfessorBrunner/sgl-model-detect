#! /usr/bin/env python
#@Author Sean McLaughlin
desc = '''
This module makes given number of "toy" images, which are actual mixtures of Gaussians with various noise filters overlaid. It also saves the true parameters
of the model in a file so they can be compared.
'''

import argparse
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('dirname', metavar = 'dirname', type = str, help = 'Directory to save the produces images to.')
parser.add_argument('N', metavar = 'N', type = int, help = 'Number of images to create.')
parser.add_argument('--noNoise', dest = 'noNoise', action = 'store_true', help =\
                                    'Save images with no noise. ')

args = parser.parse_args()

import numpy as np
from scipy.stats import poisson
from itertools import izip
from mcmcFit import mixtureOfGaussians, parseTheta, NPARAM, printTheta
import pyfits
from matplotlib import pyplot as plt
import os
import seaborn as sns
sns.set()

#The true number of Gaussians for each model
nGaussians = poisson.rvs(2, size = args.N)
nGaussians = np.where(nGaussians == 0, 1, nGaussians) #assign 0 selectiosn to 1
size = (30,30)
imageMax = 500
darkCurrent = 20 #parametrized value
chosen_cmap = 'gnuplot2'
NOISE = not args.noNoise

centerLines = []

for nImage, n in enumerate(nGaussians):
    if n == 0:
        continue

    #Make True Parameter Values
    ndim = int(n*NPARAM+2)
    trueTheta = np.zeros(ndim)
    for i in xrange(ndim):

        if i < 2: #center
            if i < n:
                j = 1
            else:
                j = 0
            mean = size[j]/2
            x = -1
            while x<0 or x>size[j]:
                x = size[j]/2+(np.random.rand()-.5)*size[j]/2

            trueTheta[i] = x

        elif i < n+2: #amp
            x = -1
            while x<0:
                x = np.random.randn()*50+imageMax/(2*(i-1))
            trueTheta[i] = x
        elif i<3*n+2: #var
            x = -1
            while x <0:
                x = np.random.randn()*6+10
            trueTheta[i] = x
        else: #corr
            x = -2
            while abs(x) > 1:
                x = np.random.randn()*.1
            trueTheta[i] = x

    yy, xx = np.indices(size)

    output = printTheta(n,trueTheta, string = True)

    #Note that c_x and c_y aren't used in this model (so 15,15 is overwritten)
    image = mixtureOfGaussians(xx,yy,15,15, trueTheta, True)

    fluxError = np.zeros(size).reshape((-1))

    for i, pixel in enumerate(image.flatten()):
        if pixel<1:
            continue
        fluxError[i] = poisson.rvs(int(pixel))
    if NOISE:
        image = fluxError.reshape(size)

    whiteNoiseMean = np.random.randn()*np.sqrt(10)+10

    if NOISE:
        image += np.random.randn(*size)*np.sqrt(whiteNoiseMean)+whiteNoiseMean #white Noise

    #Dark Current
    darkCurrentError = poisson.rvs(darkCurrent, size = size[0]*size[1]).reshape(size)-darkCurrent
    if NOISE:
        image+=darkCurrentError

    if np.any(image<=0):
        image-=image.min()*1.01 #ensure no negative or zero values
        image+=.0001

    im = plt.imshow(image, cmap = chosen_cmap)
    im.set_clim(0, image.max())
    plt.colorbar(im)
    plt.show()

    hdu = pyfits.PrimaryHDU(image)
    filename = args.dirname+'toy_image_%d'%nImage

    if 'toy_image_%d.fits'%nImage in os.listdir(args.dirname):
	    os.remove('%s.fits'%filename) #delete previous one.

    hdu.writeto('%s.fits'%filename)
    #np.savetxt(filename+'_theta', trueTheta, delimiter = ',')
    with open(filename+'_theta', 'w') as f:
        f.write(output)

    centerLines.append(' '.join(['toy_image_%d'%nImage, str(15), str(15)]))

#Gotta make a center file, too!
with open(args.dirname + 'toy_centers', 'w') as f:
    f.write('\n'.join(centerLines))
