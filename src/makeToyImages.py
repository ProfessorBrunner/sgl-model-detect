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
    lines = []
    if N == 1:
        lines.append('1 Gaussian Model\n')
    else:
        lines.append('%d Gaussians Model\n'%N)

    X, Y, As, VarXs, VarYs, Corrs = parseTheta(theta)
    for i, (a, vx, vy, p) in enumerate(izip(As, VarXs, VarYs, Corrs)):
        j = i+1
        lines.append('Gaussian %d:'%j)
        lines.append('Center %d:\t (%.2f, %.2f)'%(j,X,Y))
        lines.append('Amplitude %d:\t%.3f\nVarX %d:\t%.3f\nVarY %d:\t%.3f\nCorr %d:\t%.3f\n'%(j, a,j, vx, j, vy, j, p) )
    lines.append('\n'+'--'*20)
    return '\n'.join(lines)

np.random.seed(int(time()))

#The true number of Gaussians for each model
nGaussians = poisson.rvs(2, size = args.N)
nGaussians = np.where(nGaussians == 0, 1, nGaussians) #assign 0 selectiosn to 1
size = (30,30)
imageMax = 400
darkCurrent = 20 #parametrized value
chosen_cmap = 'gnuplot2'

centerLines = []

for nImage, n in enumerate(nGaussians):
    if n == 0:
        continue

    ndim = n*NPARAM+2
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
            #trueTheta[i] = 10**(8*np.random.rand()-4)
            #trueTheta[i] = np.random.lognormal(mean = np.log(imageMax/2), sigma = 1.0)#try logNormal near max
            x = -1
            while x<0:
                x = np.random.randn()*50+imageMax/4
            trueTheta[i] = x
            #TODO add negative in guess, or just let it explore that way?
        elif i<3*n+2: #var
            x = -1
            while x <0:
                x = np.random.randn()*6+10
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

    output= printTheta(n,trueTheta)
    print output
    c_x, c_y = trueTheta[0], trueTheta[1]
    image = sum(gaussian(xx,yy, c_x, c_y, a, varX, varY, corr) for a, varX, varY, corr in izip(*(parseTheta(trueTheta)[2:])))

    fluxError = np.zeros(size).reshape((-1))

    for i, pixel in enumerate(image.flatten()):
        if pixel<1:
            continue
        fluxError[i] = poisson.rvs(int(pixel))

    #image = fluxError.reshape(size)
    '''
    im = plt.imshow(image, cmap = chosen_cmap)
    plt.colorbar(im)
    plt.show()
    '''
    whiteNoiseMean = np.random.randn()*np.sqrt(10)+10

    #image += np.random.randn(*size)*np.sqrt(whiteNoiseMean)+whiteNoiseMean #white Noise
    '''
    im = plt.imshow(image, cmap = chosen_cmap)
    plt.colorbar(im)
    plt.show()
    '''

    #Dark Current
    darkCurrentError = poisson.rvs(darkCurrent, size = size[0]*size[1]).reshape(size)-darkCurrent
    #image+=darkCurrentError

    im = plt.imshow(image, cmap = chosen_cmap)
    im.set_clim(0, image.max())
    plt.colorbar(im)
    plt.show()


    hdu = pyfits.PrimaryHDU(image)
    filename = args.dirname+'toy_image_%d'%nImage

    if 'toy_image_%d.fits'%nImage in os.listdir(args.dirname):
	    os.remove('%s.fits'%filename)

    hdu.writeto('%s.fits'%filename)
    #np.savetxt(filename+'_theta', trueTheta, delimiter = ',')
    with open(filename+'_theta', 'w') as f:
        f.write(output)

    centerLines.append(' '.join(['toy_image_%d'%nImage, str(15), str(15)]))

#Gotta make a center file, too!
with open(args.dirname + 'toy_centers', 'w') as f:
    f.write('\n'.join(centerLines))
