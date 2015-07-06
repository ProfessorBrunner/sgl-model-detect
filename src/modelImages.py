#! /usr/bin/env python
#@Author Sean McLaughlin

desc = '''
This program models the light profile of a galaxy in a fits image, subtracts out the model and inspects the residuals for evidence of lensing.

This is the "main" module, that imports and runs the other modules in this directory
'''
import argparse
parser = argparse.ArgumentParser(description = desc)
#an option to change the number of walkers, to increase precision?
#I need a series of optional arguments to have save checkpoints along the process.
parser.add_argument('filename', metavar = 'filename', type = str, help = 'Either a fits filename or a directory of files.')
parser.add_argument('outputdir', metavar = 'output', type = str, help = 'Location to store the programs output')
parser.add_argument('imageFormat', metavar = 'format', type = str, help = 'Determines the format of images to use. Optiones are "C" for CFHTLS or "S" for SDSS.')
parser.add_argument('bands', metavar = 'bands', type = str, help = "What bands to use in the fit. If just one, will fit only that band. If 2 given, first will be\
                                                                   primary fit to and 2nd will be the one subtracted from")
parser.add_argument('centers', metavar = 'centers', type = str, help = \
    'Either a filename or a comma separate pair of coordinates for x,y. Default is to use findCenter.', default = None, nargs = '?')

parser.add_argument('--cutout', dest = 'cutout', action = 'store_const', const = True, default = False,\
                     help = 'Save a .png of the original cutout to file.')
parser.add_argument('--cutoutData', dest = 'cutoutData', action = 'store_const', const = True, default = False,\
                     help = 'Save the raw data of the cutout to file.')
parser.add_argument('--chain', dest = 'chain', action = 'store_const', const = True, default = False,\
                     help = 'Store the markov chain in the output directory. OCCUPIES A LOT OF FILESPACE.')
parser.add_argument('--subtraction', dest = 'subtraction', action = 'store_const', const = True, default = False,\
                     help = 'Store a .png of the model subtraction from the original image.')
#parser.add_argument('--residuals', dest = 'residuals', action = 'store_const', const = True, default = False,\
#                     help = 'Save a .png of the detected residuals')
parser.add_argument('--subtractionData', dest = 'subtractionData', action = 'store_const', const = True, default = False,\
                            help = 'Save the raw image data to file of the residuals.')

args = parser.parse_args()

#import matplotlib as mpl
#mpl.use('Agg')
from mcmcFit import mcmcFit, parseTheta
from nlsqFit import nlsqFit
from residualID import residualID
import imageClass
import numpy as np
from matplotlib import pyplot as plt
import seaborn
seaborn.set()
import os
from goodnessOfFit import goodnessOfFit
from itertools import izip

SHOW_IMAGES = True
chosen_cmap = 'gnuplot2'

def plotSingleImage(imageObj, bands, chosen_cmap, outputdir, name, show = False):
    for band in bands:
        fig = plt.figure()
        im = plt.imshow(imageObj.images[band], cmap = chosen_cmap)
        plt.colorbar(im)

        c_x, c_y = imageObj.center
        plt.scatter(c_x, c_y, color = 'k')

        plt.savefig(outputdir+imageObj.imageID+'_%s_%s.png'%(band, name))
        if show:
            plt.show()
        else:
            plt.clf()
            plt.close(fig)

def plotFullModel(rawImage, model, N, id, outputdir, chosen_cmap, show = False ):
    'Helper function that creates a more elaborate plot of the image and model.'
    imPlots = []
    fig = plt.figure(figsize = (30,20))
    minVal, maxVal = 0, 0

    plt.subplot(131)
    im = plt.imshow(rawImage,cmap = chosen_cmap)
    minVal = min(minVal, rawImage.min())
    maxVal = max(maxVal, rawImage.max())
    imPlots.append(im)

    plt.subplot(132)
    im = plt.imshow(model,cmap = chosen_cmap)
    minVal = min(minVal, model.min())
    maxVal = max(maxVal, model.max())
    imPlots.append(im)

    plt.subplot(133)
    calc_img = rawImage - model
    im = plt.imshow(calc_img,cmap = chosen_cmap)
    minVal = min(minVal,calc_img.min())
    maxVal = max(maxVal,calc_img.max())
    imPlots.append(im)

    for image in imPlots:
        image.set_clim(minVal, maxVal)
    cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
    fig.colorbar(imPlots[0], cax, orientation='horizontal')
    fig.suptitle('%d Gaussian Model'%N)
    plt.savefig(outputdir+id+'_%d_'%N+'fullModel.png')
    if show:
        plt.show()
    else:
        plt.clf()
        plt.close(fig)

#TODO Think of moving to mcmcFit
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

def caculateMaxGaussians(theta):
    'From the 1 G fit calculate the maximum number of Gaussians to prevent overfitting'

    #estimate how much "signal" we have, so we don't overfit
    #area enclosed within nSigma
    pixelsPerParam = 10
    nSigma = 2
    varX, varY = theta[3:5]
    #area of an ellipse
    area = np.pi*np.power(nSigma, 2)*np.sqrt(varX*varY)
    #MaxGaussians is sometimes <=2, which means only one sample is run. There should be a minimum max Gaussians.
    maxGaussians =  int(area/(4*pixelsPerParam))
    if maxGaussians < 2:
        maxGaussians = 2 #minimum value of 3
    elif maxGaussians > 10:
        maxGaussians = 10 #upper value of 10. Arbitrarily chosen and can be extended.
    print 'Max Gaussians = %d'%maxGaussians

    return maxGaussians
filename = args.filename
outputdir = args.outputdir
bands = args.bands

if len(bands) > 2 or ',' in bands:
    print 'Invalid band entry; at most 2 bands may be entered, no commas. Ex: ig'
    from sys import exit
    exit(-1)

if args.imageFormat not in ['C', 'S', 'T']:
    print 'Invalid format entry; please select "C","S" or "T"'
    from sys import exit
    exit(-1)

#determine which center finding method to use
#TODO rename useFindCenters to be more descriptive
useFindCenters = args.centers is None
if useFindCenters:
    isCoordinates = False
else:
    isCoordinates = args.centers.find(',') != -1

if isCoordinates:
    splitCenters = args.centers.split(',')
    c_x = int(splitCenters[0].strip())
    c_y = int(splitCenters[1].strip())

files = [filename, outputdir]
#Then we've been given a list of centers that will have to be read in.
if not useFindCenters and not isCoordinates:
    files.append(args.centers)
    centers = args.centers
#check files for existance

if os.path.isdir(filename) and filename[-1] != '/':
    filename+='/'

for f in files:
    if not os.path.exists(f):
        print 'ERROR: Invalied path %s'%f
        from sys import exit
        exit(-1)

primaryBand = bands[0]
secondaryBand = bands[1] if len(bands)>1 else None

#Make a dictionary of the Galaxy's image coordinates.
#TODO Rename GalaxyDict to a more descriptive name?
if not (useFindCenters or isCoordinates) :
    galaxyDict = {}
    idx = 2 if args.imageFormat == 'C' else 1 #The catalogs are formatted slightly differently
    # This was easier than fixing the catalogs, but a standard format may be preferable.
    with open(centers) as f:
        for line in f:
            splitLine = line.strip().split(' ')
            coords = [float(x) for x in splitLine[idx:]]
            galaxyDict[splitLine[0]] = coords
else:
    galaxyDict = None

imageClassDict = {'C': imageClass.CFHTLS, 'S': imageClass.SDSS, 'T': imageClass.Toy}
#the appropriate formatting for these objects
imgClass = imageClassDict[args.imageFormat]

fitter = mcmcFit#nlsqFit

imageDict = {}

#load in filenames from directory
#TODO Reads in all images, even the ones not currently being fit to. Waste of memory and time.
if os.path.isdir(filename):
    #if the input is a directory
    fileList = os.listdir(filename)
    fileDirectory = filename
    for fname in fileList:
        #Only look at .fits images
        if fname[-5:] != '.fits':
            continue
        obj = imgClass(fileDirectory+fname)

        if obj.imageID in imageDict:
            imageDict[obj.imageID].addImage(fileDirectory+fname)

        else:
            imageDict[obj.imageID] = obj

#load in lone image
#load in other bands of that image as well
else:
    obj = imgClass(filename)
    imageDict[obj.imageID] = obj
    obj.getOtherBand(bands)

for imageObj in imageDict.values():

    print 'Image ID : %s'%imageObj.imageID

    if isCoordinates:
        coords = (c_x, c_y)
    else:
        coords = None
    #adjusts the center estimate, or finds one outright if none were given.
    try:
        imageObj.calculateCenter(coords,galaxyDict)
    except KeyError:
        print 'No center for %s'%imageObj.imageID
        continue

    imageObj.cropImage()

    #Plot the Requested Cutout
    if args.cutout:
        plotSingleImage(imageObj, bands, chosen_cmap, outputdir, 'cutout', show = SHOW_IMAGES)

    if args.cutoutData:
        for band in bands:
            np.savetxt(outputdir+imageObj.imageID+'_'+band+'_cutoutData', imageObj.images[band])

    BEs = []
    prim_fits = []
    #perform first fit with 1 Gaussian
    #Then, use to charecterize max number of parameters
    c_x, c_y = imageObj.center
    print 'Fitting now'
    prim_fit,theta, be  = fitter(imageObj[primaryBand], 1, dirname = outputdir, id = imageObj.imageID, chain = args.chain)

    printTheta(1, theta)

    BEs.append(be)
    prim_fits.append(prim_fit)

    plotFullModel(imageObj[primaryBand], prim_fit, 1, imageObj.imageID, outputdir, chosen_cmap, show = SHOW_IMAGES)

    maxGaussians = caculateMaxGaussians(theta)

    #iterate until we reach our limit or BE decreases
    #maxGaussians = 4
    for N in xrange(2,maxGaussians+1):
        prim_fit, theta, be = fitter(imageObj[primaryBand], N, dirname = outputdir, id = imageObj.imageID, chain = args.chain)
        printTheta(2, theta)

        BEs.append(be)
        prim_fits.append(prim_fit)

        plotFullModel(imageObj[primaryBand],prim_fit, N, imageObj.imageID, outputdir, chosen_cmap, show = SHOW_IMAGES)

        print 'BF Diff: %.3f\t Old: %.3f\t New: %.3f'%(BEs[-1]-BEs[-2], BEs[-1], BEs[-2])
        if BEs[-1] < BEs[-2]: #new Model is worse!
        #NOTE Double-check that this is right and not supposed to be backwards
            break
            #pass

    #TODO fix scaling so it uses the calculated center rather than the image's center
    bestArg = np.argmax(np.array(BEs))
    prim_fit = prim_fits[bestArg]
    calcImgDict = {}
    calc_img = imageObj[primaryBand] - prim_fit
    calcImgDict[primaryBand] = calc_img

    c = (imageObj.center[1], imageObj.center[0])
    print 'Best model N = %d'%(bestArg+1)
    if secondaryBand is not None:
        prim_fit_scaled = prim_fit*imageObj[secondaryBand][c]/imageObj[primaryBand][c]
        calc_img = imageObj[secondaryBand] - prim_fit_scaled
        calcImgDict[secondaryBand] = calc_img

    if args.subtraction:
        plotSingleImage(imageObj, bands, chosen_cmap, outputdir, 'subtraction', show = SHOW_IMAGES)

    if args.subtractionData:
        np.savetxt(outputdir+imageObj.imageID+'_'+band+'_subtractionData', calc_img)

    goodnessOfFit(prim_fit, imageObj[primaryBand], 6*(bestArg+1), 1)

    #TODO Plotting functionality here
    #TODO Have this flagged on/off. We don't need to check for lenses on all of em.
    #TODO take advantage of calaculated center.
    #check for lens properties
    c_y, c_x = imageObj.center
    lens = residualID(calc_img,c_y, c_x )
    print 'The image %s represents a lens: %s'%(imageObj.imageID, str(lens))
