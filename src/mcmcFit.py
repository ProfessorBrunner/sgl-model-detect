#! /usr/bin/env python
#@Author Sean McLaughlin
'''
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
import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
import emcee as mc
from time import time
from matplotlib import pyplot as plt
from scipy.stats import mode, gaussian_kde
from multiprocessing import cpu_count, Pool
from itertools import izip
import seaborn
seaborn.set()


#the number of parameters per Gaussian
NPARAM = 4

def gaussian(x,y, cx, cy, a, VarX, VarY, corr):

    sigXsigY = np.sqrt(VarX*VarY)
    xDiff = x-cx
    yDiff = y-cy
    #Normally this calculation would be [xDiff, yDiff].T*np.linalg.inv(cov)*[xDiff, yDiff]
    #Multiplying through by hand yields this more efficient result.
    expVal = (-(xDiff**2)/(2*VarX)-(yDiff**2)/(2*VarY)+corr*(xDiff*yDiff)/(sigXsigY))/(1-corr**2)

    return  a*np.exp(expVal)
#TODO Make better and add more funcitons
def parseTheta(theta, movingCenter = True):
    'Helper function that splits theta into its components, for simplicity'
    if movingCenter:
        N = (len(theta)-2)/NPARAM
        X, Y = theta[0], theta[1]
        As, VarXs, VarYs, Corrs = [theta[i*N+2:(i+1)*N+2] for i in xrange(NPARAM)]
        return X, Y, As, VarXs, VarYs, Corrs
    else:
        N = len(theta)/NPARAM
        As, VarXs, VarYs, Corrs = [theta[i*N:(i+1)*N] for i in xrange(NPARAM)]
        return As, VarXs, VarYs, Corrs

def printTheta(N, theta, movingCenters = True):
    'Helper function that prints the model produced by the sampler.'
    if N == 1:
        print '1 Gaussian Model\n'
    else:
        print '%d Gaussians Model\n'%N

    if movingCenters:
        X, Y, As, VarXs, VarYs, Corrs = parseTheta(theta)
        print 'Center:\t (%.2f, %.2f)\n'%(X,Y)
    else:
        As, VarXs, VarYs, Corrs = parseTheta(theta)

    for i, (a, vx, vy, p) in enumerate(izip(As, VarXs, VarYs, Corrs)):
        j = i+1
        print 'Gaussian %d:'%j
        print 'Amplitude %d:\t%.3f\nVarX %d:\t%.3f\nVarY %d:\t%.3f\nCorr %d:\t%.3f'%(j, a,j, vx, j, vy, j, p)
        print
    print'\n'+'--'*20

#theta contains the variables for the amplitude and width
#theta = [A1,A2...An,VarX1, VarX1..., VarXN, VarY1, VarY2,...VarYN] add covs later,Cov1, Cov2,...CovN]
#if fixedCenter is false, theta has c_x, and c_y prepended.
def lnprior(theta, movingCenter, imageSize):
    #log uniform priors

    pt = parseTheta(theta, movingCenter= movingCenter)
    if movingCenter:
        N = (len(theta)-2)/NPARAM
        c_x, c_y, amps, varXs, varYs, corrs = pt
    else:
        N = len(theta)/NPARAM #save us from having to put N in the global scope
        amps, varXs, varYs, corrs = pt

    if movingCenter and any(c< 0 or c>maxC for c, maxC in izip((c_x, c_y), imageSize)):
        return -np.inf
    #TODO made negative amplitudes possible
    if any(-1e5>a or a>1e5 for a in amps):
        return -np.inf

    #enforcing order
    if any(abs(amps[i])<abs(amps[i+1]) for i in xrange(N-1)):
        return -np.inf

    if sum(amps)<0: #The whole thing has to be positive!
        return -np.inf

    #TODO find a proper bounds for this value
    for var in (varXs, varYs):
        if any(v<1e-1 or v>50 for v in var):
            return -np.inf

    if any(corr<-1 or corr>1 for corr in corrs):
        return -np.inf
    #Uniform prior
    return 1
    #log-uniform
    #lnp= -1*np.sum(np.log(theta[2*movingCenter:(NPARAM-2)*N]))
    #TODO make centerStd tunable and make an option for a uniform distrbution
    #if not movingCenter:
    #    return lnp
    #centerStd = (imageSize[0]+imageSize[1])/2 #Average size divided by 2 (~65% of the time center is within this distance of the image center)
    #return lnp - (sum(imageSize)/2-(c_x+c_y))/(2*centerStd) #logNormal for the centers


def lnlike(theta, image, xx,yy,c_x, c_y,inv_sigma2, movingCenter):

    pt = parseTheta(theta, movingCenter)

    if movingCenter:
        c_x, c_y, amps, varXs, varYs, corrs = pt
    else:
        amps, varXs, varYs, corrs = pt

    model = sum(gaussian(xx,yy,c_x,c_y,a,varX, varY, corr) for a, varX, varY, corr in izip(amps, varXs, varYs, corrs))

    diff = image-model

    #basic log normal liklihood
    #assume Gaussian errors
    #Adding Poisson weights
    return -.5*(np.sum(((diff)**2)*inv_sigma2 - np.log(inv_sigma2)))
    #return  -.5*(np.sum( (diff**2)/image-2*np.log(image)))
    #return  -.5*(np.sum( (diff**2)*image+2*np.log(image)))


#note if movingCenter is true, the c_x, c_y here are overwritten immeadiately. However, if fixedCenter is true they are needed.
def lnprob(theta, image, xx, yy, c_x, c_y, inv_sigma2, movingCenter):
    lp = lnprior(theta, movingCenter, image.shape)
    if np.isfinite(lp):
        return lp+lnlike(theta, image, xx, yy, c_x, c_y, inv_sigma2, movingCenter)
    return -np.inf

def beHelper(sample):
    #Hellper for BayesianEvidence
    #has to be top level for parallelization purposes
    logDens = np.log(beHelper.kde(sample)[0])#acquire a lock for this object?
    logp = lnlike(sample, *beHelper.args)
    return logp+np.log(beHelper.N)-logDens

#What coding standards?
beHelper.kde = None
beHelper.args = None
beHelper.N = None

def BayesianEvidence(samples, args):
    #technique taken form the code in astroML to calculate Bayesian odds. Be sure to cite
    #They use a simler method than what I employ here
    from time import time
    t0 = time()
    nCores = cpu_count()
    N,D = samples.shape

    beHelper.kde = gaussian_kde(samples.T)
    beHelper.args = args
    beHelper.N = N

    nSamples = 500
    #normalized liklihood?
    #All samples takes too long. Do a small selection
    randSamplesIdx = np.random.choice(xrange(len(samples)), size = nSamples, replace = False)
    randSamples = samples[randSamplesIdx]

    p = Pool(nCores)
    try:
        allBEs = p.map(beHelper, randSamples)
    except MemoryError:
        print 'Memory Limit Exceeded. Change the settings to better fit this setup.'
        return -np.inf
    finally:
        p.close()
        p.terminate()

    p25, p50, p75 = np.percentile(allBEs, [25, 50, 75])
    BE, dBE =  p50, 0.7413 * (p75 - p25)
    print 'BE: %f\tdBE/BE: %f'%(BE, -dBE/BE)
    #BE = np.median(allBEs)
    print 'BE Calculation time: %.6f Seconds'%(time()-t0)
    return BE

def plotChain(calc_vals,samples, n_bins, modes, means, medians, dirname, id, show = False, movingCenter = True):
    'Plot the chain of the distribution'
    fig = plt.figure(figsize = (30,15))
    fig.canvas.set_window_title('Samples %s'%id)

    ndim = len(calc_vals)

    pt = parseTheta(calc_vals, movingCenter)

    if movingCenter:
        calc_X, calc_Y, calc_as, calc_varXs, calc_varYs, calc_corrs = pt
        N = (ndim-2)/NPARAM
    else:
        calc_as, calc_varXs, calc_varYs, calc_corrs = pt
        N = ndim/NPARAM

    '''
    if N == 2: #BONUS ROUND
        chosen_idxs = np.random.choice(len(samples), size = 1e5, replace=False)
        seaborn.jointplot(samples[chosen_idxs, 1], samples[chosen_idxs,3], kind = 'kde', space = 0)#, joint_kws={'alpha': .002})
        plt.savefig(dirname + '%s_%d_scatter.png'%(id, N))
        plt.clf()
        plt.close()
    '''
    for i in xrange(ndim):
        nRows = ndim/3+(not ndim%3==0)
        plt.subplot(nRows, 3,i+1)

        offset = movingCenter*2

        if i==0 and movingCenter:
            plt.title('X : %.3f'%(calc_X))
        elif i == 1 and movingCenter:
            plt.title('Y : %.3f'%(calc_Y))
        elif i < N+offset:
            j = i-offset
            plt.title("Amplitude %d: %.3f"%(j+1, calc_as[j]))
        elif i < 2*N+offset:
            j = i-N-offset
            plt.title("$\sigma^2_X$ %d: %.3f"%(j+1, calc_varXs[j]))
        elif i < 3*N+offset:
            j = i-2*N-offset
            plt.title("$\sigma^2_Y$ %d: %.3f"%(j+1, calc_varYs[j]))
        else:
            j = i-3*N-offset
            plt.title('$\\rho$ %d: %.3f'%(j+1, calc_corrs[j]))

        plt.hist(samples[:, i], bins = n_bins)
        plt.vlines(modes[i],0,5e3, colors = ['r'], label = 'mode')
        plt.vlines(means[i],0,5e3, colors = ['g'], label = 'mean')
        plt.vlines(medians[i],0,5e3, colors = ['m'],label = 'median')
        plt.legend()

    if dirname is not None and id is not None:
        plt.savefig(dirname + '%s_%d_chain.png'%(id,N))
    if show:
        plt.show()
    else:
        plt.clf()
        plt.close(fig)

#Centers should still be needed for initial guess
#TODO change the order of these parameters around, it doesn't make intuitive sense anymore. Make sure to change nlsq, too!
def mcmcFit(image, N, c_x, c_y, movingCenters, n_walkers = 2000, dirname = None, id = None, chain = False):
    np.random.seed(int(time()))
    #TODO Check if i'm going to exceed memory limits?
    t0 = time()
    #numpy arrays of the indicies, used in the calculations
    yy, xx = np.indices(image.shape)

    #error used in the liklihood. Its value does not seem to change the results much.
    #TODO Will change resutls depending on iamge scale!
    #Represents the std of the error, which is assumed Gaussian
    inv_sigma2 = pow(10, 3)

    #parameters for the emcee sampler.
    ndim = N*NPARAM #1 Amplitude and 3 Radial dimentions
    if movingCenters:
        ndim+=2
    nsteps = 400 #Sample more for larger models
    nburn = int(nsteps*.25)

    #initial guess
    pos = []
    imageMax = image.max()

    for walk in xrange(n_walkers):
        row = np.zeros(ndim)
        if movingCenters:
            for i in xrange(2):
                mean = image.shape[i]/2
                x = -1
                while x<0 or x>image.shape[i]:
                    #x = mean*np.random.randn()+mean
                    x = np.random.rand()*image.shape[i]
                row[i] = x

        for i in xrange(2*movingCenters,ndim):
            if i < 2*movingCenters+N: #amp
                #row[i] = 10**(8*np.random.rand()-4)
                row[i] = np.random.lognormal(mean = np.log(imageMax/2), sigma = 1.0)#try logNormal near max
                #TODO add negative in guess, or just let it explore that way?
            elif i<ndim-N: #var
                #TODO fix not fitting other Gaussians
                row[i] = 10**(2*np.random.rand()-1)
                #row[i] = 10**(8*np.random.rand()-4)
            else: #corr
                #row[i] = 2*np.random.rand()-1
                #Trying a normal dist. rather and a uniform.
                #The assumption being Galaxies are more likely to be spherical than very elliptical
                x = -2
                while abs(x) > 1:
                    x = np.random.randn()
                row[i] = x
        pos.append(row)

    args = (image, xx, yy, c_x, c_y, inv_sigma2, movingCenters)
    sampler = mc.EnsembleSampler(n_walkers, ndim, lnprob, args = args,threads = cpu_count())
    #run the sampler. Longest running line in the code
    sampler.run_mcmc(pos, nsteps)

    samples = sampler.chain[:,nburn:,:].reshape((-1, ndim))
    sampler.pool.terminate()#there's a bug in emcee that creates daemon threads. This kills them.
    del(sampler)
    #save samples to file
    if chain and dirname is not None and id is not None:
        np.savetxt(dirname+'_%s_%d_chain'%(id,N), samples, delimiter = ',')

    #TODO MAP as chain median?
    #the MAP is simply the mean of the chain
    calc_means = samples.mean(axis = 0)

    p25, calc_medians, p75 = np.percentile(samples,[25,50,75], axis = 0)
    spread = .7413 * (p75 - p25)
    #for i in xrange(ndim):
    #   print calc_vals[i], spread[i]

    n_bins = int(np.sqrt(samples.shape[0]))#arbitrary

    #calculate the mode from a histogram
    calc_modes = np.zeros(ndim)
    #calc_vals = np.zeros(ndim)

    for i in xrange(ndim):
        #if (movingCenters and 1< i < ndim-N) or i<ndim-N:
        #    hist, bin_edges = np.histogram(np.log10(samples[:,i]), bins = n_bins)
        #    max_idx = np.argmax(hist)
        #   #NOTE unsure if I should take the mean over the exponents or values
        #    calc_modes[i] = 10**((bin_edges[max_idx]+bin_edges[max_idx+1])/2)#center of peak
        #else:
        hist, bin_edges = np.histogram(samples[:,i], bins = n_bins)
        max_idx = np.argmax(hist)
        calc_modes[i] = (bin_edges[max_idx]+bin_edges[max_idx+1])/2

    calc_vals = calc_medians
    #calc_modes = calc_vals

    plotChain(calc_vals,samples, n_bins, calc_modes, calc_means, calc_medians, dirname, id, show = True, movingCenter = movingCenters)

    if movingCenters:
        cx, cy, calc_as, calc_varXs, calc_varYs, calc_corrs = parseTheta(calc_vals, movingCenter= movingCenters)
    else:
        calc_as, calc_varXs, calc_varYs, calc_corrs = parseTheta(calc_vals, movingCenter= movingCenters)

    calc_img = sum(gaussian(xx,yy,cx,cy,a,varX, varY, corr) for a, varX, varY, corr in izip(calc_as, calc_varXs, calc_varYs, calc_corrs))

    #Calculate the evidence for this model
    print 'Fitting Time: %.4f Seconds'%(time()-t0)
    print 'Chain Size: %d Megabytes'%(samples.nbytes/(1024**2))
    BE = BayesianEvidence(samples, args)
    return calc_img, calc_vals, BE

