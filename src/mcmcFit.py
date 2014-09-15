#! /usr/bin/env python
#@Author Sean McLaughlin
desc ='''
This module fits to an image using MCMC techniques, specifically using the package emcee. 
This module contains the function fitImage. It can be run as main or imported.

mcmcFit(image, N, c_x, c_y, n_walkers = 600, ddof = 0)
image: numpy array of the image
N: the number of Gaussians to use in the fit
c_x, c_y: the center of the object in the image.
n_walkers: the number of walkers, default 600
ddof: the change in degrees of freedom, for the chiSquare test

return:
calc_img: the image with the residual subtracted out
chi2stat: the chi2 statistic
p: the probability of the model. Usually 0.

To run as main:
python mcmcFit.py filename 2 100 100

with optional n_walkers and ddof
'''
import numpy as np
import emcee as mc
from time import time
from matplotlib import pyplot as plt
from scipy.stats import mode,chisquare
from multiprocessing import cpu_count
from sklearn.cluster import KMeans

def gaussian(x,y, cx, cy, a, cov):
    pos = np.array([x,y])
    mu = np.array([cx,cy])
    muMPos = pos-mu
    invCov = np.linalg.inv(cov)
    return a*np.exp(np.dot(muMPos.T, np.dot(invCov, muMPos)))
    
#theta contains the variables for the amplitude and width
#theta = [A1,A2...An,R11, R12, R13, R14, R21,...Rn4]
def lnprior(theta):
    #log uniform priors
    #simplified down to uniform priors, but these represent the exponents not the values themselves
    N = len(theta)/5 #save us from having to put N in the global scope
    amps = theta[:N]
    rads = theta[N:]
    #NOTE to solve the uniqueness problem I can require the amplitdues are in order. Slower, but all will converge to separe values.
    if not all(-1<a<3 for a in amps):
        return -np.inf
    #My scheme for constructing a matrix requires they be uniform
    if not all(0<r<1 for r in rads):
        return -np.inf
    return 0

def lnlike(theta, image, xx,yy,c_x, c_y,inv_sigma2):
    N = len(theta)/5
    amps = theta[:N]
    rads = theta[N:]
    rads = np.array(rads).reshape(N, 4)
    covs = []
    for rad in rads:
        #construct a positive definite matrix
        #TODO make into a function cuz it'll be copy pasted a lot.
        cov = 0.5 - rad.reshape((2, 2))
        cov = np.triu(cov)
        cov += cov.T - np.diag(cov.diagonal())
        covs.append( np.dot(cov,cov))

    model = np.zeros(image.shape) 
    for a,r in zip(amps, covs):
        if a<0: #if a<0 makes the amplitdue 0.
            continue 
        model+=gaussian(xx,yy,c_x, c_y, 10**a, cov)

    diff = image-model
    #optioal punishment for overfitting vs. underfitting
    #largerThanImage = diff<0
    #diff[largerThanImage]*=10 #punish going over more than under

    #basic log normal liklihood
    return -.5*(np.sum(((diff)**2)*inv_sigma2-np.log(inv_sigma2)))

def lnprob(theta, image, xx, yy, c_x, c_y, inv_sigma2):
    lp = lnprior(theta)
    if np.isfinite(lp):
        return lp+lnlike(theta, image, xx, yy, c_x, c_y, inv_sigma2)
    return -np.inf

def calcsCluster(samples, N, decimals = 2, n_clusters = 3):
    #select the parameters with a clustering approach

    #FIX for covariance matricies
    sampled_as = samples[:, :N].reshape((-1))
    sampled_rs = samples[:, N:].reshape((-1))

    data = np.c_[sampled_as, sampled_rs]
    #reshape the data such that it's 2-D, with amps on one axis and radii on the other.

    n_clusters = N+1
    #fit to clusters
    k_means = KMeans(init = 'k-means++',n_clusters = n_clusters, n_init = 50)
    labels = k_means.fit_predict(data)

    #round the data for binning and mode selection
    roundData = np.round_(data, decimals = decimals)
    clusters = [roundData[labels == i] for i in xrange(n_clusters)]

    #get a list of the highest points and their average count
    allModes = []
    for cluster in clusters:
           point, counts = mode(cluster)
           allModes.append((point[0][0], point[0][1], counts.mean()))

    allModes.sort(key = lambda x:x[-1], reverse = True) #sort by highest count
    allModes = np.array(allModes)
    #select the 1st N popular points in the clusters
    modes = allModes[:N,:2]

    calc_as, calc_rs = modes[:,0], modes[:,1]
    return calc_as, calc_rs 

def mcmcFit(image, N, c_x, c_y, n_walkers = 600, ddof = 0, filename = None):

    np.random.seed(int(time()))

    img_y, img_x = image.shape

    #numpy arrays of the indicies, used in the calculations
    yy, xx = np.indices(image.shape)

    #error used in the liklihood. It's value does not seem to change the results much.
    err = .1
    inv_sigma2 = 1./(err**2)

    #parameters for the emcee sampler.
    ndim = N*5 #1 Amp and 5 Rad dimentions
    nburn = int(n_walkers*.1)
    nsteps = 200

    #initial guess
    pos = []
    for walk in xrange(n_walkers):
        row = np.zeros(ndim)
        for n in xrange(N):
            row[n] = 4*np.random.rand()-1
            for i in xrange(4):
                row[n+N+i] = np.random.rand()
        pos.append(row)

    #sometimes the center is not exactly accurate. This part finds the maximum in the region around the center.
    dy, dx = np.unravel_index(image[c_y-1:c_y+2, c_x-1:c_x+2].argmax(), (3,3))
    dy,dx = dy-1, dx-1
    c_y, c_x = c_y+dy, c_x+dx

    args = (image, xx, yy, c_x, c_y, inv_sigma2) 
    sampler = mc.EnsembleSampler(n_walkers, ndim, lnprob, args = args, \
                                                    threads = cpu_count()) 
    #run the sampler. Longest line in the code
    sampler.run_mcmc(pos, nsteps)
    samples = sampler.chain[:,nburn:,:].reshape((-1, ndim))
    sampler.pool.terminate()#there's a bug in emcee that creates daemon threads. This kills them.
    del(sampler)
    #save samples to file
    if filename is not None:
        np.savetxt(filename, samples, delimiter = ',')

    #using clustering on the samples to select the parameters from the posterior
    #TODO fix toward covariance
    calc_as, calc_rs = calcsCluster(samples, N)

    calc_img = sum(gaussian(xx,yy,c_x,c_y,10**a,10**r) for a,r in zip(calc_as, calc_rs))
    #calcuate the chi2 test
    ddof = -2*N + ddof 
    chi2stat, p = chisquare(image, f_exp = calc_img, ddof = ddof,axis = None)

    return calc_img, chi2stat, p

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
    parser.add_argument('ddof', metavar = 'ddof', type = int, help = 'Change in the degree of freedom',nargs = '?', default = 0)

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
    calc_img, chi2stat, p = mcmcFit(image, N, c_x, c_y, args.n_walkers, args.ddof) 
    plt.subplot(121)
    im = plt.imshow(image)
    plt.colorbar()
    plt.subplot(122)
    im = plt.imshow(image-calc_img)
    plt.colorbar()
    plt.show()
