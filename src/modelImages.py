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
    parser.add_argument('imageFormat', metavar = 'format', type = str, help = 'Determines the format of images to use. Optiones are "C" for CFHTLS or "S" for SDSS.')
    parser.add_argument('bands', metavar = 'bands', type = str, help = "What bands to use in the fit. If just one, will fit only that band. If 2 given, first will be\
                                                                       primary fit to and 2nd will be the one subtracted from")
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

    bands = args.bands

    if args.imageFormat not in ['C', 'S']:
        print 'Invalid format entry; please select "C" or "S"'
        from sys import exit
        exit(-1)

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

    if os.path.isdir(filename) and filename[-1] != '/':
        filename+='/'
    
    for f in files:
        if not os.path.exists(f):
            print 'ERROR: Invalied path %s'%f
            from sys import exit
            exit(-1)
    #begin setting up the input dict
    #This is a straightforward way to carry around various options
    #TODO Reconsider the input dict; it's become more nightmare than tool
    inputDict = {}
    inputDict['filename'] = filename
    inputDict['output']= outputdir
    inputDict['useFindCenters'] = useFindCenters
    inputDict['isCoordinates'] = isCoordinates
    inputDict['isDir'] = os.path.isdir(filename)
    inputDict['primaryBand'] = bands[0]
    inputDict['secondaryBand'] = bands[1] if len(bands)>1 else None

    if isCoordinates:
        inputDict['coords'] = (cx, cy)

    #Make a dictionary of the Galaxy's image coordinates.
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
    from mcmcFit import mcmcFit
    from residualID import residualID
    import imageClass
    import numpy as np

    imageClassDict = {'C': imageClass.CFHTLS, 'S': imageClass.SDSS}
    #the appropriate formatting for these objects
    imgClass = imageClassDict[args.imageFormat]
    imageDict = {}

    #load in filenames from directory
    if inputDict['isDir']:
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
    else:
        obj = imgClass(filename)
        imageDict[obj.imageID] = obj 

    for imageObj in imageDict.itervalues():
        #savefile name for sample chain
        name = inputDict['output']+baseName+'_samples' if inputDict['chain'] else None
        coords = None
        galaxyDict = None

        if inputDict['isCoordinates']:
            coords = inputDict['coords']

        if 'galaxyDict' in inputDict:
            galaxyDict = inputDict['galaxyDict']

        imageObj.calculateCenters(coords,galaxyDict)
        #print 'Cutout'
        imageObj.cropImage(inputDict['cutout'], inputDict['output'])

        BFs = []
        prim_fits = []
        #perform first fit with 1 Gaussian
        #Then, use to charecterize max number of parameters
        c_x, c_y = imageObj.center
        prim_fit,theta, bf  = mcmcFit(imageObj[inputDict['primaryBand']], 1, c_x, c_y, filename = name)
        BFs.append(bf)
        prim_fits.append(prim_fit)

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
            prim_fit, theta, bf = mcmcFit(imageObj[inputDict['primaryBand']], n, c_x, c_y, filename = name)
            BFs.append(bf)
            prim_fits.append(prim_fit)
            if BFs[-1]/BFs[-1] > 1: #new Model is worse!
                break

        bestArg = np.argmax(np.array(BFs))
        print 'Best model N = %d'%(bestArg+1)
        if inputDict['secondaryBand'] is not None:
            prim_fit = prim_fits[bestArg]
            prim_fit_scaled = prim_fit*imageObj[inputDict['secondaryBand']][c]/imageObj[inputDict['primaryBand']][c]
            calc_img = imageObj[inputDict['secondaryBand']] - prim_fit_scaled

        else:
            calc_img = imageObj[inputDict['primaryBand']] - prim_fit
        '''
        from matplotlib import pyplot as plt
        print 'Model'
        im = plt.imshow(prim_fit)
        plt.colorbar(im)
        plt.show()
        plt.clf()
        plt.close()
        '''
        if inputDict['subtraction']:
            from matplotlib import pyplot as plt
            im = plt.imshow(calc_img)
            plt.colorbar(im)
            plt.savefig(inputDict['output']+imageObj.imageID+'_subtraction.png')
            #print 'Subtraction'
            #plt.show()
            plt.clf()
            plt.close()

        if inputDict['subtractionData']:
            import numpy as np
            np.savetxt(inputDict['output']+imageObj.imageID+'_residualData', calc_img)

        #TODO Plotting functionality here
        #check for lens properties
        lens = residualID(calc_img, c_x, c_y)
        print 'The image %s represents a lens: %s'%(imageObj.imageID, str(lens))

if __name__ == '__main__':
    main()
