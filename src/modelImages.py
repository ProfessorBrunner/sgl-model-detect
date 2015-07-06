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

parser.add_argument('--fixedCenters', dest = 'fixedCenters', action = 'store_const', const = True, default = False,\
                     help = 'Use the given centers as fixed values. Otherwise, they will be treated as initial guesses')
#TODO should this be an input argument instead of what I've got here?
parser.add_argument('--nlsq', dest = 'nlsq', action = 'store_const', const = True, default = False,\
                    help = 'Instead of the standard MCMC, use a NLSQ method to fit a MOG to the image.')
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
from mcmcFit import mcmcFit, parseTheta, printTheta
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

def plotSingleImage(imageDict, bands, chosen_cmap, outputdir, name, show = False):

    for band in bands:
        fig = plt.figure()
        plImg = imageDict[band]
        im = plt.imshow(plImg, cmap = chosen_cmap)
        plt.colorbar(im)

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

def caculateMaxGaussians(theta, movingCenters = True):
    'From the 1 G fit calculate the maximum number of Gaussians to prevent overfitting'

    #estimate how much "signal" we have, so we don't overfit
    #area enclosed within nSigma
    pixelsPerParam = 10
    nSigma = 2
    varX, varY = theta[1+2*movingCenters:3+2*movingCenters]
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

fitter = mcmcFit#nlsqFit#mcmcFit

if args.nlsq:
    fitter = nlsqFit

imageClassDict = {'C': imageClass.CFHTLS, 'S': imageClass.SDSS, 'T': imageClass.Toy}
#the appropriate formatting for these objects
imgClass = imageClassDict[args.imageFormat]
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
    #savefile name for sample chain
    print '-'*30
    print '*'*30
    print '_'*30
    print 'Image ID : %s'%imageObj.imageID

    #make a folder for all outputs
    #All outputs for this object will go in the same folder!
    imOutputDir = outputdir+imageObj.imageID+'/'
    if os.path.isdir(imOutputDir):
        for f in os.listdir(imOutputDir):
            os.remove(imOutputDir+f)
    else:
        os.mkdir(imOutputDir)

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
        plotSingleImage(imageObj.images, bands, chosen_cmap, imOutputDir, 'cutout', show = SHOW_IMAGES)

    if args.cutoutData:
        for band in bands:
            np.savetxt(imOutputDir+imageObj.imageID+'_'+band+'_cutoutData', imageObj.images[band])

    stats = []
    prim_fits = []
    #perform first fit with 1 Gaussian
    #Then, use to charecterize max number of parameters
    c_x, c_y = imageObj.center
    print 'Fitting now'
    prim_fit,theta, stat  = fitter(imageObj[primaryBand], 1, c_x, c_y, not args.fixedCenters, dirname = imOutputDir, id = imageObj.imageID, chain = args.chain)
    printTheta(1, theta, movingCenters= not args.fixedCenters)

    stats.append(stat)
    prim_fits.append(prim_fit)

    plotFullModel(imageObj[primaryBand], prim_fit, 1, imageObj.imageID, imOutputDir, chosen_cmap, show = SHOW_IMAGES)

    maxGaussians = caculateMaxGaussians(theta)

    #TODO delete
    #maxGaussians = 4
    for n in xrange(2,maxGaussians+1):
        prim_fit, theta, stat = fitter(imageObj[primaryBand], n, c_x, c_y,not args.fixedCenters, dirname = imOutputDir, id = imageObj.imageID, chain = args.chain)

        printTheta(2, theta, movingCenters= not args.fixedCenters)

        if np.isnan(stat):
            print 'NaN raised for Gaussian %d'%n
            break

        stats.append(stat)
        prim_fits.append(prim_fit)

        plotFullModel(imageObj[primaryBand], prim_fit, n, imageObj.imageID, imOutputDir, chosen_cmap, show = SHOW_IMAGES)

        print ' Old: %.3f\t New: %.3f'%(stats[-2], stats[-1])


        if args.nlsq:
            #do a chi2 test
            test = abs(stats[-1]-1)> abs(stats[-2]-1)
        else:
            test = stats[-1]<stats[-2] #equivalent to BF < 1

        if test: #new Model is worse!
        #NOTE Double-check that this is right and not supposed to be backwards
            break
            #pass

    #TODO fix scaling so it uses the calculated center rather than the image's center
    if args.nlsq:
        bestArg = np.argmin(np.abs(np.array(stats)-1))
    else:
        bestArg = np.argmax(np.array(stats))

    prim_fit = prim_fits[bestArg]
    calcImgDict = {}
    calc_img = imageObj[primaryBand] - prim_fit
    calcImgDict[primaryBand] = calc_img
    print 'Best model N = %d'%(bestArg+1)
    if secondaryBand is not None:
        c = (imageObj.center[1], imageObj.center[0])
        prim_fit_scaled = prim_fit*imageObj[secondaryBand][c]/imageObj[primaryBand][c]
        calc_img = imageObj[secondaryBand] - prim_fit_scaled
        calcImgDict[secondaryBand] = calc_img

    if args.subtraction:
        plotSingleImage(calcImgDict, bands, chosen_cmap, imOutputDir, 'subtraction', show = SHOW_IMAGES)

    if args.subtractionData:
        np.savetxt(imOutputDir+imageObj.imageID+'_'+band+'_residualData', calc_img)

    #TODO Plotting functionality here
    #TODO Have this flagged on/off. We don't need to check for lenses on all of em.
    #check for lens properties
    goodnessOfFit(prim_fit, imageObj[primaryBand], 4*(bestArg+1)+2*(not args.fixedCenters), 1)
    #TODO take advantage of calculated center
    c_y, c_x = imageObj.center
    #lens = residualID(calc_img, c_x, c_y)
    #print 'The image %s represents a lens: %s'%(imageObj.imageID, str(lens))
