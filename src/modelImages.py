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
    #TODO Adjust input for center finding. Centers can be given as hard centers or initial guesses. Perhaps a flag?
    #Some guess has to be entered because cropping won't be possible otherewise
    parser.add_argument('centers', metavar = 'centers', type = str, help = \
        'Either a filename or a comma separate pair of coordinates for x,y. Default is to use findCenter.', default = None, nargs = '?')

    parser.add_argument('--fixedCenters', dest = 'fixedCenters', action = 'store_const', const = True, default = False,\
                         help = 'Use the given centers as fixed values. Otherwise, they will be treated as initial guesses')
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

    if len(bands) > 2 or ',' in bands:
        print 'Invalid band entry; at most 2 bands may be entered, no commas. Ex: ig'
        from sys import exit
        exit(-1)

    if args.imageFormat not in ['C', 'S']:
        print 'Invalid format entry; please select "C" or "S"'
        from sys import exit
        exit(-1)

    #determine which center finding method to use
    #TODO Is findCenters still useful at all? Probably no sense deleting it.
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
    #Then we've been given a list of centers that will have to be read in.
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

    inputDict['fixedCenters'] = args.fixedCenters
    inputDict['cutout'] = args.cutout
    inputDict['cutoutData'] = args.cutoutData
    inputDict['chain'] = args.chain
    inputDict['subtraction'] = args.subtraction
    #inputDict['residuals']=args.residuals
    inputDict['subtractionData'] = args.subtractionData

    if isCoordinates:
        inputDict['coords'] = (cx, cy)

    #Make a dictionary of the Galaxy's image coordinates.
    #TODO Rename GalaxyDict to a more descriptive name?
    elif not useFindCenters:
        galaxyDict = {}
        idx = 2 if args.imageFormat == 'C' else 1 #The catalogs are formatted slightly differently
        # This was easier than fixing the catalogs, but a standard format may be preferable.
        with open(centers) as f:
            for line in f:
                splitLine = line.strip().split(' ')
                coords = [float(x) for x in splitLine[idx:]]
                galaxyDict[splitLine[0]] = coords

        inputDict['galaxyDict'] = galaxyDict


    from mcmcFit import mcmcFit
    from residualID import residualID
    import imageClass
    import numpy as np

    imageClassDict = {'C': imageClass.CFHTLS, 'S': imageClass.SDSS}
    #the appropriate formatting for these objects
    imgClass = imageClassDict[args.imageFormat]
    imageDict = {}

    #load in filenames from directory
    #TODO Reads in all images, even the ones not currently being fit to. Waste of memory and time.
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
    #load in other bands of that image as well
    else:
        obj = imgClass(filename)
        imageDict[obj.imageID] = obj
        obj.getOtherBand(bands)

    #For now, randomize so I see different values while testing.
    #TODO Delete after testing, so it's more stable.
    x = imageDict.values()
    np.random.shuffle(x)
    for imageObj in x:
        #savefile name for sample chain
        name = inputDict['output']+imageObj.imageID+'_samples' if inputDict['chain'] else None
        print 'Image ID : %s'%imageObj.imageID
        coords = None
        galaxyDict = None

        if inputDict['isCoordinates']:
            coords = inputDict['coords']

        if 'galaxyDict' in inputDict:
            galaxyDict = inputDict['galaxyDict']

        #adjusts the center estimate, or finds one outright if none were given.
        imageObj.calculateCenter(coords,galaxyDict)

        imageObj.cropImage()

        #Plot the Requested Cutout
        #TODO Where is cutout data?
        if inputDict['cutout']:
            for band in bands:
                from matplotlib import pyplot as plt
                plt.figure()
                im = plt.imshow(imageObj.images[band])
                plt.colorbar(im)
                c_x, c_y = imageObj.center
                plt.scatter(c_x, c_y, color = 'k')
                #If leaving at this spot, fix filename
                #plt.savefig(filename)
                plt.show()
                plt.clf()
                plt.close()

        if inputDict['cutoutData']:
            import numpy as np
            for band in bands:
                np.savetxt(inputDict['output']+imageObj.imageID+'_cutoutData', imageObj.images[band])

        BEs = []
        prim_fits = []
        #perform first fit with 1 Gaussian
        #Then, use to charecterize max number of parameters
        c_x, c_y = imageObj.center
        print 'Fitting now'
        prim_fit,theta, be  = mcmcFit(imageObj[inputDict['primaryBand']], 1, c_x, c_y, not inputDict['fixedCenters'], filename = name)
        print 'Gaussian #1'
        print theta
        print'--'*15
        BEs.append(be)
        prim_fits.append(prim_fit)

        #Necessary because imageObj.center is in the opposite order of numpy slicing.
        if inputDict['fixedCenters']:
            c = (int(c_y), int(c_x))
        else:
            c = theta[1], theta[0]
        print c

        #TODO Delete, only for testing

        prim_fit_scaled = prim_fit*imageObj[inputDict['secondaryBand']][c]/imageObj[inputDict['primaryBand']][c]
        calc_img = imageObj[inputDict['secondaryBand']] - prim_fit_scaled
        from matplotlib import pyplot as plt
        im = plt.imshow(calc_img)
        plt.colorbar(im)
        plt.title('%d'%1)
        plt.show()

        #estimate how much "signal" we have, so we don't overfit
        #area enclosed within nSigma
        pixelsPerParam = 10
        nSigma = 2

        varX, varY = theta[1:3]
        #area of an ellipse
        area = np.pi*np.power(nSigma, 2)*np.sqrt(varX*varY)
        #MaxGaussians is sometimes <=2, which means only one sample is run. There should be a minimum max Gaussians.
        maxGaussians =  int(area/(4*pixelsPerParam))
        if maxGaussians < 2:
            maxGaussians = 2 #minimum value of 3
        elif maxGaussians > 10:
            maxGaussians = 10 #upper value of 10. Arbitrarily chosen and can be extended.
        print 'Max Gaussians = %d'%maxGaussians

        #iterate until we reach our limit or BE decreases
        for n in xrange(2,maxGaussians+1):
            prim_fit, theta, be = mcmcFit(imageObj[inputDict['primaryBand']], n, c_x, c_y,not inputDict['fixedCenters'], filename = name)
            print 'Gaussian #%d'%n
            print theta
            print '--'*15

            BEs.append(be)
            prim_fits.append(prim_fit)
            #TODO Delete, only for testing

            prim_fit_scaled = prim_fit*imageObj[inputDict['secondaryBand']][c]/imageObj[inputDict['primaryBand']][c]
            calc_img = imageObj[inputDict['secondaryBand']] - prim_fit_scaled

            im = plt.imshow(calc_img)
            plt.colorbar(im)
            plt.title('%d'%n)
            plt.show()

            print 'Diff: %.3f\t Old: %.3f\t New: %.3f'%(BEs[-1]-BEs[-2], BEs[-1], BEs[-2])
            if BEs[-1] < BEs[-2]: #new Model is worse!
            #NOTE Double-check that this is right and not supposed to be backwards
                break

        #TODO fix scaling so it uses the calculated center rather than the image's center
        bestArg = np.argmax(np.array(BEs))
        print 'Best model N = %d'%(bestArg+1)
        if inputDict['secondaryBand'] is not None:
            prim_fit = prim_fits[bestArg]
            prim_fit_scaled = prim_fit*imageObj[inputDict['secondaryBand']][c]/imageObj[inputDict['primaryBand']][c]
            calc_img = imageObj[inputDict['secondaryBand']] - prim_fit_scaled

        else:
            calc_img = imageObj[inputDict['primaryBand']] - prim_fit

        #TODO Delete; only for testing
        from matplotlib import pyplot as plt
        print 'Model'
        im = plt.imshow(prim_fit)
        plt.title('Model, N = %d'%(bestArg+1))
        plt.colorbar(im)
        plt.show()
        plt.clf()
        plt.close()

        plt.title('BEs')
        plt.plot(BEs)
        plt.scatter(bestArg, BEs[bestArg], color = 'r')
        plt.show()
        plt.clf()
        plt.close()

        if inputDict['subtraction']:
            from matplotlib import pyplot as plt
            im = plt.imshow(calc_img)
            plt.colorbar(im)
            plt.savefig(inputDict['output']+imageObj.imageID+'_subtraction.png')
            print 'Subtraction'
            plt.show()
            plt.clf()
            plt.close()

        if inputDict['subtractionData']:
            import numpy as np
            np.savetxt(inputDict['output']+imageObj.imageID+'_residualData', calc_img)

        #TODO Plotting functionality here
        #TODO Have this flagged on/off. We don't need to check for lenses on all of em.
        #check for lens properties
        lens = residualID(calc_img, c[1], c[0])
        print 'The image %s represents a lens: %s'%(imageObj.imageID, str(lens))

if __name__ == '__main__':
    main()
