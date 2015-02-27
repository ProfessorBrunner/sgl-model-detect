#! /usr/bin/env python
#@Author Sean McLaughlin

def main():
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
    parser.add_argument('centers', metavar = 'centers', type = str, help = 'Either a filename or a comma separate pair of coordinates for x,y. Default is to use findCenter.', default = None, nargs = '?')
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
    filename = args.filename
    outputdir = args.outputdir

    #determine which center finding method to use
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
    #check files for existance
    import os

    if os.isdir(filename) and filename[-1] != '/':
        filename+='/'
    
    for f in files:
        if not os.path.exists(f):
            print 'ERROR: Invalied path %s'%f
            from sys import exit
            exit(-1)
    #being setting up the input dict
    #This is a straightforward way to carry around various options
    inputDict = {}
    inputDict['filename'] = filename
    inputDict['output']= outputdir
    inputDict['useFindCenters'] = useFindCenters
    inputDict['isCoordinates'] = isCoordinates
    inputDict['isDir'] = os.path.isdir(filename)

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
    inputDict['subtraction'] = args.subtraction
    #inputDict['residuals']=args.residuals
    inputDict['subtractionData'] = args.subtractionData

    from cropImage import cropImage
    from findCenter import findCenter
    from mcmcFit import mcmcFit
    from residualID import residualID
    import numpy as np
    import pyfits

    #load in filenames from directory
    if inputDict['isDir']:
        fileList = os.listdir(filename)
        baseNames = set()
        trackedObjs = []
        for fname in fileList:
            if fname[:6] == 'CFHTLS':
                baseNames.add(fname[:-7])

        baseNames = list(baseNames)
        fileDirectory = filename
    #load in lone image
    else:
        baseName = filename[:-7]
        lineIndex = baseName.rfind('/')
        fileDirectory, baseName = baseName[:lineIndex], baseName[lineIndex:]
        baseNames = [baseName]

    bands = ['g', 'i']
    for baseName in baseNames:
        #savefile name for sample chain
        name = inputDict['output']+baseName+'_samples' if inputDict['chain'] else None

        images = {}
        for band in bands:
            #load in image and find its center through designated method
            fitsImage = pyfits.open(fileDirectory+baseName+'_'+band+'.fits')
            image = fitsImage[0].data
            if inputDict['useFindCenters']:
                c_y, c_x = findCenter(image)
            elif inputDict['isCoordinates']:
                c_x, c_y = inputDict['coords']
            else:
                c_x, c_y = inputDict['galaxyDict'][baseName[7:]]
            #crop down to size, and possibly save
            image, c_x, c_y = cropImage(image, c_x, c_y, plot = inputDict['cutout'], filename = inputDict['output'] + baseName+'_'+band+'_cutout.png')
            if inputDict['cutoutData']:
                import numpy as np
                np.savetxt(inputDict['output']+baseName+'_'+band+'_cutoutData', image)

            images[band] = image

        BFs = []
        i_fits = []
        #perform first fit with 1 Gaussian
        #Then, use to charecterize max number of parameters
        i_fit,theta, bf  = mcmcFit(images['i'], 1, c_x, c_y, filename = name)
        BFs.append(bf)
        i_fits.append(i_fit)

        c = (int(c_y), int(c_x))

        #estimate how much "signal" we have, so we don't overfit
        #area enclosed within nSigma
        pixelsPerParam = 10
        nSigma = 2

        varX, varY = theta[1:3]
        #area of an ellipse
        area = np.pi*np.power(nSigma, 2)*np.sqrt(varX*varY)

        maxGaussians =  int(area/(4*pixelsPerParam))

        #iterate until we reach our limit or BF decreases
        for n in xrange(2,maxGaussians):
            i_fit, theta, bf = mcmcFit(images['i'], n, c_x, c_y, filename = name)
            BFs.append(bf)
            i_fits.append(i_fit)
            if BFs[-1]/BFs[-1] > 1: #new Model is worse!
                break

        bestArg = np.argmax(np.array(BFs))
        print 'Best model N = %d'%(bestArg+1)
        i_fit = i_fits[bestArg]
        i_fit_scaled = i_fit*images['g'][c]/images['i'][c]
        calc_img = images['g'] - i_fit_scaled

        if inputDict['subtraction']:
            from matplotlib import pyplot as plt
            im = plt.imshow(calc_img)
            plt.colorbar(im)
            plt.savefig(inputDict['output']+baseName+'_subtraction.png')
            plt.show()
            plt.clf()
            plt.close()

        if inputDict['subtractionData']:
            import numpy as np
            np.savetxt(inputDict['output']+baseName+'_residualData', calc_img)

        #TODO Plotting functionality here
        #check for lens properties
        lens = residualID(calc_img, c_x, c_y)
        print 'The image %s represents a lens: %s'%(baseName, str(lens))

if __name__ == '__main__':
    main()
