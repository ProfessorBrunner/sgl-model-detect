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
    '''
    minVal, maxVal = 0,0
    from matplotlib import pyplot as plt
    imList = []
    plt.subplot(221)
    im = plt.imshow(data, cmap = 'gnuplot2')
    minVal = min(minVal, data.min())
    maxVal = max(maxVal, data.max())
    imList.append(im)
    plt.subplot(222)
    im = plt.imshow(1/np.sqrt(model), cmap = 'gnuplot2')
    minVal = min(minVal, model.min())
    maxVal = max(maxVal, model.max())
    imList.append(im)
    plt.subplot(223)
    calc_img = data-model
    im = plt.imshow(calc_img, cmap = 'gnuplot2')
    minVal = min(minVal, calc_img.min())
    maxVal = max(maxVal, calc_img.max())
    imList.append(im)
    plt.subplot(224)
    im = plt.imshow(errImage, cmap = 'gnuplot2')
    #minVal = min(minVal, errImage.min())
    #maxVal = max(maxVal, errImage.max())
    imList.append(im)

    for image in imList:
        image.set_clim(minVal, maxVal)

    plt.colorbar(imList[1])

    plt.show()
    '''

    chi2Statistic = np.sum( i**2 for i in errImage.reshape((-1)))
    chi2dof = chi2Statistic/dof
    print 'chi2dof = %.3f'%chi2dof
    return chi2dof

    #rv = chi2(dof)
    #x = np.linspace(0, 2*dof, num = 500)
    #plt.plot(x, rv.pdf(x))
    #plt.scatter(chi2Statistic, rv.pdf(chi2Statistic), color = 'r')
    #plt.show()