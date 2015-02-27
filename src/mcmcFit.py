#! /usr/bin/env python
#@Author Sean McLaughlin
desc ='''
This module fits to an image using MCMC techniques, specifically using the package emcee. 
This module contains the function fitImage. It can be run as main or imported.

mcmcFit(image, N, c_x, c_y, n_walkers = 600, filename = None)
image: numpy array of the image
N: the number of Gaussians to use in the fit
c_x, c_y: the center of the object in the image.
n_walkers: the number of walkers, default 600
filename: If not None, the location to store the sample chain

return:
calc_img: the image with the residual subtracted out
calc_vals: The parameter vector corresponding to the image
BF: The Bayes Factor for this model

To run as main:
python mcmcFit.py filename 2 100 100

with optional n_walkers and saveFile
'''
import numpy as np
import emcee as mc
from time import time
from matplotlib import pyplot as plt
from scipy.stats import mode, gaussian_kde
from multiprocessing import cpu_count
from itertools import izip

#TODO: there's a problem with the vectorization of this calculation. Isn't efficient in this form.

#the number of parameters per Gaussian
NPARAM = 4
def gaussian(x,y, cx, cy, a, cov):
    muMPos = np.dstack([x-cx, y-cy]) 
    invCov = np.linalg.inv(cov)
    output = np.zeros(x.shape)
    for i, row in enumerate(muMPos):
        for j, vec in enumerate(row):
            output[i,j] =np.dot(vec, np.dot(invCov, vec))  
        
    return a*np.exp(-.5*output) 

#theta contains the variables for the amplitude and width
#theta = [A1,A2...An,VarX1, VarX1..., VarXN, VarY1, VarY2,...VarYN] add covs later,Cov1, Cov2,...CovN]
#TODO Put number of parameters in the global scope
def lnprior(theta):
    #log uniform priors
    N = len(theta)/NPARAM #save us from having to put N in the global scope
    amps, varXs, varYs, corrs = [theta[i*N:(i+1)*N] for i in xrange(NPARAM)]

    if any(1e-1>a or a>1e3 for a in amps):
        return -np.inf

    #enforcing order
    if any(amps[i]<amps[i+1] for i in xrange(N-1)):
        return -np.inf

    #TODO find a proper bounds for this value
    for var in (varXs, varYs):
        if any(v<1e-2 or v>1e3 for v in var):
            return -np.inf

    if any(corr<-1 or corr>1 for corr in corrs):
        return -np.inf

    #log Uniform prior
    return -1*np.sum(np.log(theta[:(NPARAM-1)*N]))

def lnlike(theta, image, xx,yy,c_x, c_y,inv_sigma2):
    N = len(theta)/NPARAM
    amps, varXs, varYs, corrs = [theta[i*N:(i+1)*N] for i in xrange(NPARAM)]

    covariance_mats = []
    for varX, varY, corr in izip(varXs, varYs, corrs):
        #construct a covariance matrix
        cov = corr*np.sqrt(varX*varY)
        mat = np.array([varX, corr, corr, varY]).reshape((2,2))
        covariance_mats.append(mat)

    model = np.zeros(image.shape) 
    for a,cov in izip(amps, covariance_mats):
        model+=gaussian(xx,yy,c_x, c_y, a, cov)

    diff = image-model

    #basic log normal liklihood
    #assume Gaussian errors
    return -.5*(np.sum(((diff)**2)*inv_sigma2-np.log(inv_sigma2)))

def lnprob(theta, image, xx, yy, c_x, c_y, inv_sigma2):
    lp = lnprior(theta)
    if np.isfinite(lp):
        return lp+lnlike(theta, image, xx, yy, c_x, c_y, inv_sigma2)
    return -np.inf

def BayesFactor(samples, theta, args):
    #technique taken form the code in astroML to calculate Bayesian odds. Be sure to cite
    #They use a simler method than what I employ here
    N,D = samples.shape

    kde = gaussian_kde(samples.T)
    logDens = kde(theta)[0]
    logp = lnlike(theta, args)

    BF = logp+np.log(N)-logDens

    return BF

def mcmcFit(image, N, c_x, c_y, n_walkers = 1000, filename = None):

    np.random.seed(int(time()))

    #numpy arrays of the indicies, used in the calculations
    yy, xx = np.indices(image.shape)

    #error used in the liklihood. It's value does not seem to change the results much.
    #Represents the std of the error, which is assumed Gaussian
    inv_sigma2 = pow(.1, -2)

    #parameters for the emcee sampler.
    ndim = N*NPARAM #1 Amp and 3 Rad dimentions
    nsteps = 200
    nburn = int(nsteps*.25)

    #initial guess
    pos = []
    for walk in xrange(n_walkers):
        row = np.zeros(ndim)
        for i in xrange(ndim):
            if i < N: #amp
                row[i] = 10**(4*np.random.rand()-1)
            elif i<ndim-N: #var
                row[i] = 10**(5*np.random.rand()-2)
            else: #corr
                row[i] = 2*np.random.rand()-1
        pos.append(row)
    #TODO Try removing this section and see what happens
    #sometimes the center is not exactly accurate. This part finds the maximum in the region around the center.
    dy, dx = np.unravel_index(image[c_y-1:c_y+2, c_x-1:c_x+2].argmax(), (3,3))
    dy,dx = dy-1, dx-1
    c_y, c_x = c_y+dy, c_x+dx

    args = (image, xx, yy, c_x, c_y, inv_sigma2) 
    sampler = mc.EnsembleSampler(n_walkers, ndim,\
                                 lnprob, args = args,threads = cpu_count())
    #run the sampler. Longest line in the code
    sampler.run_mcmc(pos, nsteps)

    samples = sampler.chain[:,nburn:,:].reshape((-1, ndim))
    sampler.pool.terminate()#there's a bug in emcee that creates daemon threads. This kills them.
    del(sampler)
    #save samples to file
    if filename is not None:
        np.savetxt(filename, samples, delimiter = ',')

    #the MAP is simply the mean of the chain
    #calc_vals = samples.mean(axis = 0)

    n_bins = int(np.sqrt(samples.shape[0])/4)#arbitrary
    calc_vals = np.zeros(ndim)

    for i in xrange(ndim):
        if i< ndim-N:
            hist, bin_edges = np.histogram(np.log10(samples[:,i]), bins = n_bins)
            max_idx = np.argmax(hist)
            #NOTE unsure if I should take the mean over the exponents or values
            calc_vals[i] = 10**((bin_edges[max_idx]+bin_edges[max_idx+1])/2)#center of peak
        else:
            hist, bin_edges = np.histogram(samples[:,i], bins = n_bins)
            max_idx = np.argmax(hist)
            calc_vals[i] = (bin_edges[max_idx]+bin_edges[max_idx+1])/2

    calc_as, calc_varXs, calc_varYs, calc_corrs = [calc_vals[i*N:(i+1)*N] for i in xrange(NPARAM)]

    covariance_mats = []
    for varX, varY, corr in izip(calc_varXs, calc_varYs, calc_corrs):
        #construct a covariance matrix
        cov = corr*np.sqrt(varX*varY)
        mat = np.array([varX, cov, cov, varY]).reshape((2,2))
        covariance_mats.append(mat)

    calc_img = sum(gaussian(xx,yy,c_x,c_y,a,cov) for a,cov in izip(calc_as, covariance_mats))
    BF = BayesFactor(samples, calc_vals, args)
    return calc_img, calc_vals, BF

if __name__ == '__main__':
    import argparse
    import pyfits
    from cropImage import cropImage 

    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('filename', metavar = 'fname', type = str, help = 'The name of the fits file to be read in.')
    parser.add_argument('nGaussians', metavar = 'N', type = int, help = 'Number of Gaussians to use in the fit.', default = 2)
    parser.add_argument('center_x', metavar = 'cx', type = int, help = 'The center in x')
    parser.add_argument('center_y', metavar = 'cy', type = int, help = 'The center in y')
    parser.add_argument('n_walkers', metavar = 'n_walkers', type = int, help = 'Number of walkers',nargs = '?', default = 1000)
    parser.add_argument('saveFile', metavar = 'sfile', type = str, help = 'Where to store the sample chain.',nargs = '?', default = 0)

    args = parser.parse_args()

    filename = args.filename

    try:
        fitsImage = pyfits.open(filename)
    except IOError:
        print 'ERROR: Invalid filename.'
        from sys import exit
        exit(-1)

    image = fitsImage[0].data
    #TODO Make findcenter find the center of the images if none are passed in.
    c_y, c_x = args.center_y, args.center_x

    image, c_x, c_y = cropImage(image, c_x, c_y)

    N = args.nGaussians
    calc_img, theta, BF = mcmcFit(image, N, c_x, c_y, args.n_walkers, args.saveFile)
    plt.subplot(121)
    im = plt.imshow(image)
    plt.colorbar()
    plt.subplot(122)
    im = plt.imshow(image-calc_img)
    plt.colorbar()
    plt.show()
