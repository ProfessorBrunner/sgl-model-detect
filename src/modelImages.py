#! /usr/bin/env python
#@Author Sean McLaughlin
desc = '''
This program models the light profile of a galaxy in a fits image, subtracts out the model and inspects the residuals for evidence of lensing.

This is the "main" module, that imports and runs the other modules in this directory
'''
import argparse
parser = argparse.ArgumentParser(description = desc)
#the number of Gaussians to use? (I think I may want to use a chi2 test to find the best one.)
#an option to change the number of walkers, to increase precision?
#I need a series of optional arguments to have save checkpoints along the process.
#Among those I need:
#Store the plot of the ID'd residuals
parser.add_argument('filename', metavar = 'filename', type = str, help = 'Either a fits filename or a directory of files.')
parser.add_argument('outputdir', metavar = 'output', type = str, help = 'Location to store the programs output')
parser.add_argument('centers', metavar = 'centers', type = str, help = 'Either a filename or a comma separate pair of coordinates for x,y. Default is to use findCenter.', default = None, nargs = '?')
parser.add_argument('--cutout', dest = 'cutout', action = 'store_const', const = True, default = False,\
                     help = 'Save a .png of the original cutout to file.')
parser.add_argument('--cutoutData', dest = 'cutoutData', action = 'store_const', const = True, default = False,\
                     help = 'Save the raw data of the cutout to file.')
parser.add_argument('--chain', dest = 'chain', action = 'store_const', const = True, default = False,\
                     help = 'Store the markov chain in the output directory. OCCUPIES A LOT OF FILESPACE.')
parser.add_argument('--triangle', dest = 'triangle', action = 'store_const', const = True, default = False,\
                     help = 'Save a .png of a triangle plot of the chain.')
parser.add_argument('--subtraction', dest = 'subtraction', action = 'store_const', const = True, default = False,\
                     help = 'Store a .png of the model subtraction from the original image.')
parser.add_argument('--residuals', dest = 'residuals', action = 'store_const', const = True, default = False,\
                     help = 'Save a .png of the detected residuals')
parser.add_argument('--residualData', dest = 'residualData', action = 'store_const', const = True, default = False,\
                            help = 'Save the raw image data to file of the residuals.')
parser.add_argument('--bands', dest = 'bands', action = 'store_const', const = True, default = False,\
                    help = 'For one file, assume that an image of a similar name but of a different band is in the same directory. Fits to i but subtracts from g.')

args = parser.parse_args()
import os
filename = args.filename
outputdir = args.outputdir

useFindCenters = args.centers is None
if useFindCenters:
    isCoordinates = False
else:
    isCoordinates = args.centers.find(',') != -1

if isCoordinates:
    splitCenters = args.centers.split(',')
    cx = int(splitCenters[0].strip())
    cy = int(splitCenters[1].strip())

files = [filename, outputdir]
if not useFindCenters and not isCoordinates:
    files.append(args.centers)
    centers = args.centers

for f in files:
    if not os.path.exists(f):
        print 'ERROR: Invalied path %s'%f
        from sys import exit
        exit(-1)

inputDict = {'filename':filename, 'output': outputdir,'useFindCenters':useFindCenters, 'isCoordinates':isCoordinates, 'isDir':os.path.isdir(filename)}

if isCoordinates:
    inputDict['coords'] = (cx, cy)

elif not useFindCenters:
    galaxyDict = {}
    with open(centers) as f:
        for line in f:
            splitLine = line.strip().split(' ')
            coords = [float(x) for x in splitLine[2:]]
            galaxyDict[splitLine[0]] = coords

    inputDict['galaxyDict'] = galaxyDict

inputDict['cutout'] = args.cutout
inputDict['cutoutData'] = args.cutoutData
inputDict['chain'] = args.chain
inputDict['triangle'] = args.triangle
inputDict['subtraction'] = args.subtraction
inputDict['residuals']=args.residuals
inputDict['residualData'] = args.residualData
inputDict['bands'] = args.bands

from cropImage import cropImage
from findCenter import findCenter
from mcmcFit import mcmcFit
from residualID import residualID
import pyfits

if not inputDict['isDir']:
#its a file, fit on one image
    baseName = filename[:-7] 
    lineIndex = baseName.rfind('/')
    fileDirectory, baseName = baseName[:lineIndex], baseName[lineIndex:]
    name = inputDict['output']+baseName+'_samples' if inputDict['chain'] else None
    if inputDict['bands']: #do multiple bands
        bands = ['g', 'i']
        images = {}
        for band in bands:
            fitsImage = pyfits.open(fileDirectory+baseName+'_'+band+'.fits')
            image = fitsImage[0].data
            if inputDict['useFindCenters']:
                c_y, c_x = findCenter(image)
            elif inputDict['isCoordinates']:
                c_x, c_y = inputDict['coords']
            else:
                c_x, c_y = inputDict['galaxyDict'][baseName[8:]]
            print 'Original centers: (%d, %d)'%(c_x, c_y)
            image, c_x, c_y = cropImage(image, c_x, c_y, plot = inputDict['cutout'], filename = inputDict['output'] + baseName+'_'+band+'_cutout.png')
            if inputDict['cutoutData']:
                import numpy as np
                np.savetxt(inputDict['output']+baseName+'_'+band+'_cutoutData', image)
            images[band] = image    
        #TODO Fix ddof so chi2stat is correct!
        print 'Cropped centers: (%d, %d)'%(c_x, c_y)
        i_fit, i_stat, i_p = mcmcFit(images['i'], 2, c_x, c_y, filename = name)
        c = (int(c_y), int(c_x))
        i_fit = i_fit*images['g'][c]/images['i'][c]
        calc_img = images['g'] - i_fit
    else:
        fitsImage = pyfits.open(filename)
        image = fitsImage[0].data

        if inputDict['useFindCenters']:
            c_y, c_x = findCenter(image)
        elif inputDict['isCoordinates']:
            c_x, c_y = inputDict['coords']
        else:
            c_x, c_y = inputDict['galaxyDict'][baseName[8:]]
        image, c_x, c_y = cropImage(image, c_x, c_y, plot = inputDict['cutout'], filename = inputDict['output']+baseName+'_cutout.png')
        if inputDict['cutoutData']:
            import numpy as np
            np.savetxt(inputDict['output']+baseName+'_cutoutData', image)
        #TODO Fix ddof so chi2stat is correct!
        img_fit, chi2stat, p = mcmcFit(image,3, c_x, c_y, filename = name)
        calc_img = image-img_fit
    if inputDict['residualData']:
        import numpy as np
        np.savetxt(inputDict['output']+baseName+'_residualData', calc_img)

    if inputDict['subtraction']:
        from matplotlib import pyplot as plt
        im = plt.imshow(calc_img)
        plt.colorbar(im)
        plt.scatter(c_x, c_y, color = 'm')
        plt.show()
        #plt.savefig(inputDict['output']+baseName+'_subtraction.png')
        plt.close()

    #TODO Plotting functionality here
    lens = residualID(image, c_x, c_y)
    print lens

#a directory is handled differenly than a single file
else :
    fileList = os.listdir(filename)
    baseNames = set()
    trackedObjs = []
    for fname in fileList:
        if fname[:6] == 'CFHTLS':
            baseNames.add(fname[:-7])

    baseNames = list(baseNames)
    for baseName in baseNames:
        lineIndex = baseName.rfind('/')
        fileDirectory, baseName = baseName[:lineIndex], baseName[lineIndex:]
        bans = ['g', 'i']
        for band in bands:
            fitsImage = pyfits.open(fileDirectory + basename + '_'+ band + '.fits')
            fullImage =fitsImage[0].data
            image, c_x, 
        name = inputDict['output']+baseName+'_samples' if inputDict['chain'] else None
        bands = ['g', 'i']
        images = {}
        for band in bands:
            fitsImage = pyfits.open(filename+baseName+'_'+band+'.fits')
            image = fitsImage[0].data
            if inputDict['useFindCenters']:
                c_y, c_x = findCenter(image)
            elif inputDict['isCoordinates']:
                c_x, c_y = inputDict['coords']
            else:
                c_x, c_y = inputDict['galaxyDict'][baseName[8:]]

            image, c_x, c_y = cropImage(image, c_x, c_y, plot = inputDict['cutout'], filename = inputDict['output'] + baseName+'_'+band+'_cutout.png')
            if inputDict['cutoutData']:
                import numpy as np
                np.savetxt(inputDict['output']+baseName+'_'+band+'_cutoutData', image)
            images[band] = image    
        #TODO Fix ddof so chi2stat is correct!
        i_fit, i_stat, i_p = mcmcFit(images['i'], 3, c_x, c_y, filename = name)
        c = (int(c_y), int(c_x))
        i_fit = i_fit*images['g'][c]/images['i'][c]
        calc_img = images['g'] - i_fit
        if inputDict['residualData']:
            import numpy as np
            np.savetxt(inputDict['output']+baseName+'_residualData', calc_img)

        if inputDict['subtraction']:
            print 'Here'
            from matplotlib import pyplot as plt
            im = plt.imshow(calc_img)
            plt.colorbar(im)
            plt.show()
            #plt.savefig(inputDict['output']+baseName+'_subtraction.png')
            plt.close()

        #TODO Plotting functionality here
        lens = residualID(image, c_x, c_y)
        print lens

