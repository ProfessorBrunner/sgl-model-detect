#! /usr/bin/env python
#@Author Sean McLaughlin
'''
This model has the same API as mcmcFit.py. It takes the same inputs, and returns the same values (except Be will always be 1, because evidence doens't have meaning
in this context). It performs a non-linear least squares fit to the data. It will be used to compare different fitting techniques.
'''

from scipy.optimize import curve_fit
from mcmcFit import gaussian, parseTheta
import numpy as np
from goodnessOfFit import goodnessOfFit
from itertools import izip

NPARAM = 4

def f(x,*args):
    x = x.reshape((30,30,2))
    xx, yy = x[:,:,0], x[:,:,1]

    c_x, c_y = args[0], args[1]
    calc_img = sum(gaussian(xx,yy, c_x, c_y, a, varX, varY, corr) for a, varX, varY, corr in izip(*(parseTheta(args)[2:])))

    return calc_img.flatten()

def nlsqFit(image, N, c_x, c_y, movingCenters, n_walkers = 2000, dirname = None, id = None, chain = False):
    #note that n_walkers and chain will not be used at all by this function. Just to comply with the other API

    yy, xx = np.indices(image.shape)

    ydata = image.flatten() #curve_fit takes only flat arrays. The functions I write will reshape.
    xdata = np.dstack((xx,yy)).flatten()

    ndim = N*NPARAM #1 Amplitude and 3 Radial dimentions
    if movingCenters:
        ndim+=2

    #initial guess
    pos = np.zeros(ndim)
    imageMax = image.max()

    if movingCenters:
        for i in xrange(2):
            mean = image.shape[i]/2
            x = -1
            while x<0 or x>image.shape[i]:
                #x = mean*np.random.randn()+mean
                x = np.random.rand()*image.shape[i]
            pos[i] = x

    for i in xrange(2*movingCenters,ndim):
        if i < 2*movingCenters+N: #amp
            #pos[i] = 10**(8*np.random.rand()-4)
            pos[i] = np.random.lognormal(mean = np.log(imageMax/2), sigma = 1.0)#try logNormal near max
            print pos[i]
            #TODO add negative in guess, or just let it explore that way?
        elif i<ndim-N: #var
            #TODO fix not fitting other Gaussians
            pos[i] = 10**(2*np.random.rand()-1)
            #pos[i] = 10**(8*np.random.rand()-4)
        else: #corr
            #pos[i] = 2*np.random.rand()-1
            #Trying a normal dist. rather and a uniform.
            #The assumption being Galaxies are more likely to be spherical than very elliptical
            x = -2
            while abs(x) > 1:
                x = np.random.randn()
            pos[i] = x

    popt, pcov = curve_fit(f,xdata, ydata, p0 = pos, maxfev = int(1e6) )

    c_x, c_y = popt[0], popt[1]
    calc_img = sum(gaussian(xx,yy, c_x, c_y, a, varX, varY, corr) for a, varX, varY, corr in izip(*(parseTheta(popt)[2:])))

    chi2 = goodnessOfFit(calc_img, image, ndim, pow(10,-2))

    return calc_img, popt, chi2