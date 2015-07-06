#! /usr/bin/env python
#@Author Sean McLaughlin
'''
This model has the same API as mcmcFit.py. It takes the same inputs, and returns the same values (except Be will always be 1, because evidence doens't have meaning
in this context). It performs a non-linear least squares fit to the data. It will be used to compare different fitting techniques.
'''

from scipy.optimize import curve_fit
from mcmcFit import gaussian, parseTheta, printTheta
import numpy as np
from goodnessOfFit import goodnessOfFit
from itertools import izip

NPARAM = 4

imageShape = (30,30)

def f(x,*args):
    global imageShape
    x = x.reshape((imageShape[0],imageShape[1],2))
    xx, yy = x[:,:,0], x[:,:,1]

    c_x, c_y = args[0], args[1]
    calc_img = sum(gaussian(xx,yy, c_x, c_y, a, varX, varY, corr) for a, varX, varY, corr in izip(*(parseTheta(args)[2:])))

    return calc_img.flatten()

def nlsqFit(image, N, c_x, c_y, movingCenters, n_walkers = 2000, dirname = None, id = None, chain = False):
    #note that n_walkers and chain will not be used at all by this function. Just to comply with the other API

    yy, xx = np.indices(image.shape)
    global imageShape
    imageShape = image.shape

    ydata = image.flatten() #curve_fit takes only flat arrays. The functions I write will reshape.
    xdata = np.dstack((xx,yy)).flatten()

    ndim = N*NPARAM #1 Amplitude and 3 Radial dimentions
    if movingCenters:
        ndim+=2

    #initial guess
    pos = np.zeros(ndim)
    imageMax = image.max()

    #Try deterministic initial guess
    if movingCenters:
        pos[0] = 15
        pos[1] = 15
    for i in xrange(2*movingCenters, ndim):
        if i<2*movingCenters+N:#amp
            pos[i] = imageMax/N
        elif i<ndim-N: #var
            j = i-2*movingCenters-N if i<2*movingCenters+2*N else i-2*movingCenters-2*N
            pos[i] = 15/((j+1)*N)+1
        else: #corr
            pos[i] = 0

    try:
        popt, pcov = curve_fit(f,xdata, ydata, p0 = pos, maxfev = int(1e5) )
    except RuntimeError: #exceeded maxfev
        print 'Max Evals exceeded for %d Gaussians'%N
        return np.zeros(image.shape), np.zeros(ndim), np.NaN

    c_x, c_y = popt[0], popt[1]
    calc_img = sum(gaussian(xx,yy, c_x, c_y, a, varX, varY, corr) for a, varX, varY, corr in izip(*(parseTheta(popt)[2:])))

    chi2 = goodnessOfFit(calc_img, image, ndim, pow(10,-2))

    return calc_img, popt, chi2