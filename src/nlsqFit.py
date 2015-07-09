#! /usr/bin/env python
#@Author Sean McLaughlin
'''
This model has the same API as mcmcFit.py. It takes the same inputs, and returns the same values (except it returns a Chi2 stat, because evidence doens't have meaning
in this context). It performs a non-linear least squares fit to the data. It will be used to compare different fitting techniques.
'''

from scipy.optimize import curve_fit
from mcmcFit import mixtureOfGaussians
import numpy as np

NPARAM = 4

#Global variables so they can be accessed by the f function without being passed in as arguments.
imageShape = (30,30)
c = (15,15)
_movingCenter = True

def goodnessOfFit(model, data,k, sigma):
    SSR = np.sum(p**2 for p in (data-model).flatten())
    print 'SSR = %e'%SSR #Sum of Squared Residuals
    errImage = (data-model)/np.sqrt(model) #Poisson noise, sigma = sqrt(mu)
    N =1
    for dim in errImage.shape:
        N*=dim
    #N data points in the image
    dof = N - k

    chi2Statistic = np.sum( i**2 for i in errImage.reshape((-1)))
    chi2dof = chi2Statistic/dof
    print 'chi2dof = %.3f'%chi2dof
    return chi2dof

def f(x,*args):
    global imageShape
    x = x.reshape((imageShape[0],imageShape[1],2))
    xx, yy = x[:,:,0], x[:,:,1]

    c_x, c_y = c[1], c[0]
    
    model = mixtureOfGaussians(xx, yy, c_x, c_y, args, movingCenter= _movingCenter )

    return model.flatten()

def nlsqFit(image, N, c_x, c_y, movingCenter, n_walkers = 2000, dirname = None, id = None, chain = False):
    #note that n_walkers and chain will not be used at all by this function. Just to comply with the other API

    yy, xx = np.indices(image.shape)
    global imageShape
    imageShape = image.shape

    if not movingCenter:
        global _movingCenter
        _movingCenter = False
        global c
        c = (c_y, c_x)

    ydata = image.flatten() #curve_fit takes only flat arrays. The functions I write will reshape.
    xdata = np.dstack((xx,yy)).flatten()

    ndim = N*NPARAM #1 Amplitude and 3 Radial dimentions
    if movingCenter:
        ndim+=2

    #initial guess
    pos = np.zeros(ndim)
    imageMax = image.max()

    #Deterministic initial guess
    if movingCenter:
        pos[0] = 15
        pos[1] = 15
    for i in xrange(2*movingCenter, ndim):
        if i<2*movingCenter+N:#amp
            pos[i] = imageMax/N
        elif i<ndim-N: #var
            j = i-2*movingCenter-N if i<2*movingCenter+2*N else i-2*movingCenter-2*N
            pos[i] = 15/((j+1)*N)+1
        else: #corr
            pos[i] = 0

    try:
        popt, pcov = curve_fit(f,xdata, ydata, p0 = pos, maxfev = int(1e5) )
    except RuntimeError: #exceeded maxfev
        print 'Max Evals exceeded for %d Gaussians'%N
        return np.zeros(image.shape), np.zeros(ndim), np.NaN


    model = mixtureOfGaussians(xx, yy, c_x, c_y, popt, movingCenter)

    chi2 = goodnessOfFit(model, image, ndim, pow(10,-2))

    return model, popt, chi2