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

args = parser.parse_args()

import numpy as np
from scipy.stats import poisson, uniform
from itertools import izip
from mcmcFit import gaussian, parseTheta, NPARAM
import pyfits
from matplotlib import pyplot as plt
from time import time
import os
import seaborn as sns
sns.set()

def printTheta(N, theta):
    'Helper function that prints the model produced by the sampler.'
    if N == 1:
        print '1 Gaussian Model\n'
    else:
        print '%d Gaussians Model\n'%N

    Xs, Ys, As, VarXs, VarYs, Corrs = parseTheta(theta)
    for i, (x, y, a, vx, vy, p) in enumerate(izip(Xs, Ys, As, VarXs, VarYs, Corrs)):
        j = i+1
        print 'Gaussian %d:'%j
        print 'Center %d:\t (%.2f, %.2f)'%(j,x,y)
        print 'Amplitude %d:\t%.3f\nVarX %d:\t%.3f\nVarY %d:\t%.3f\nCorr %d:\t%.3f'%(j, a,j, vx, j, vy, j, p)
        print
    print'\n'+'--'*20

np.random.seed(int(time()))

#The true number of Gaussians for each model
nGaussians = poisson.rvs(2, size = args.N)
size = (30,30)
imageMax = 400
darkCurrent = 20 #parametrized value
chosen_cmap = 'gnuplot2'

centerLines = []

for nImage, n in enumerate(nGaussians):
    if n == 0:
        continue

    ndim = n*NPARAM
    trueTheta = np.zeros(ndim)
    for i in xrange(ndim):

        if i < 2*n: #center
            if i < n:
                j = 1
            else:
                j = 0
            mean = size[j]/2
            x = -1
            while x<0 or x>size[j]:
                x = size[j]/2+(np.random.rand()-.5)*size[j]/2
            trueTheta[i] = x

        elif i < 3*n: #amp
            #trueTheta[i] = 10**(8*np.random.rand()-4)
            #trueTheta[i] = np.random.lognormal(mean = np.log(imageMax/2), sigma = 1.0)#try logNormal near max
            x = -1
            while x<0:
                x = np.random.randn()*100+imageMax/2
            trueTheta[i] = x
            #TODO add negative in guess, or just let it explore that way?
        elif i<5*n: #var
            x = -1
            while x <0:
                x = np.random.randn()*10+15
            trueTheta[i] = x
            #trueTheta[i] = 10**(8*np.random.rand()-4)
        else: #corr
            #trueTheta[i] = 2*np.random.rand()-1
            #Trying a normal dist. rather and a uniform.
            #The assumption being Galaxies are more likely to be spherical than very elliptical
            x = -2
            while abs(x) > 1:
                x = np.random.randn()*.1
            trueTheta[i] = x

    yy, xx = np.indices(size)

    printTheta(n,trueTheta)

    image = sum(gaussian(xx,yy,c_x, c_y, a, varX, varY, corr) for c_x, c_y, a, varX, varY, corr in izip(*parseTheta(trueTheta)))

    fluxError = np.zeros(size).reshape((-1))

    for i, pixel in enumerate(image.flatten()):
        if pixel<1:
            continue
        fluxError[i] = poisson.rvs(int(pixel))

    image = fluxError.reshape(size)
    '''
    im = plt.imshow(image, cmap = chosen_cmap)
    plt.colorbar(im)
    plt.show()
    '''
    whiteNoiseMean = np.random.randn()*np.sqrt(30)+30
    print whiteNoiseMean

    image += np.random.randn(*size)*np.sqrt(whiteNoiseMean)+whiteNoiseMean #white Noise
    '''
    im = plt.imshow(image, cmap = chosen_cmap)
    plt.colorbar(im)
    plt.show()
    '''

    #Dark Current
    darkCurrentError = poisson.rvs(darkCurrent, size = size[0]*size[1]).reshape(size)-darkCurrent
    image+=darkCurrentError

    im = plt.imshow(image, cmap = chosen_cmap)
    plt.colorbar(im)
    plt.show()


    hdu = pyfits.PrimaryHDU(image)
    filename = args.dirname+'toy_image_%d'%nImage

    if 'toy_image_%d.fits'%nImage in os.listdir(args.dirname):
	    os.remove('%s.fits'%filename)

    hdu.writeto('%s.fits'%filename)
    np.savetxt(filename+'_theta', trueTheta, delimiter = ',')

    centerLines.append(' '.join(['toy_image_%d'%nImage, str(15), str(15)]))

#Gotta make a center file, too!
with open(args.dirname + 'toy_centers', 'w') as f:
    f.write('\n'.join(centerLines))
