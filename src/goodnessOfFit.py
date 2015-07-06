#! /usr/bin/env python
#@Author Sean Mclaughlin

desc ='''
This file contains function(s) for testing the goodness of the MCMC fit. As opposed to residualID, this program tests that
the modeling of an LRG worked like it was supposed to, and doesn't look for anything other than that.

'''
from scipy.stats import chi2
import numpy as np
from matplotlib import pyplot as plt

def goodnessOfFit(model, data,k, sigma):
    errImage = (data-model)/np.sqrt(model) #Poisson noise, sigma = sqrt(mu)
    N =1
    for dim in errImage.shape:
        N*=dim
    #N data points in the image
    dof = N - k

    chi2Statistic = np.sum( i**2 for i in errImage.reshape((-1)))
    chi2dof = chi2Statistic/dof
    print chi2dof
    return chi2dof

    rv = chi2(dof)
    x = np.linspace(0, 2*dof, num = 500)
    #plt.plot(x, rv.pdf(x))
    #plt.scatter(chi2Statistic, rv.pdf(chi2Statistic), color = 'r')
    #plt.show()