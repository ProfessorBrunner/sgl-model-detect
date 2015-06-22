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
    import os

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

    from mcmcFit import mcmcFit
    from residualID import residualID
    import imageClass
    import numpy as np
    import seaborn
    seaborn.set()

    imageClassDict = {'C': imageClass.CFHTLS, 'S': imageClass.SDSS}
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

    #For now, randomize so I see different values while testing.
    #TODO Delete after testing, so it's more stable.
    x = imageDict.values()
    np.random.shuffle(x)
    chosen_cmap = 'jet'
    for imageObj in x:
        #savefile name for sample chain
        name = outputdir+imageObj.imageID+'_samples' if args.chain else None
        print 'Image ID : %s'%imageObj.imageID

        if isCoordinates:
            coords = (c_x, c_y)
        else:
            coords = None
        #adjusts the center estimate, or finds one outright if none were given.
        imageObj.calculateCenter(coords,galaxyDict)

        imageObj.cropImage()

        #Plot the Requested Cutout
        if args.cutout:
            for band in bands:
                from matplotlib import pyplot as plt
                plt.figure()
                im = plt.imshow(imageObj.images[band], cmap = chosen_cmap)
                plt.colorbar(im)
                c_x, c_y = imageObj.center
                plt.scatter(c_x, c_y, color = 'k')
                plt.savefig(outputdir+imageObj.imageID+'_'+band+'_cutout.png')
                plt.show()
                plt.clf()
                plt.close()

        if args.cutoutData:
            import numpy as np
            for band in bands:
                np.savetxt(outputdir+imageObj.imageID+'_'+band+'_cutoutData', imageObj.images[band])

        BEs = []
        prim_fits = []
        #perform first fit with 1 Gaussian
        #Then, use to charecterize max number of parameters
        c_x, c_y = imageObj.center
        print 'Fitting now'
        prim_fit,theta, be  = mcmcFit(imageObj[primaryBand], 1, c_x, c_y, not args.fixedCenters, filename = name)
        print 'Gaussian #1'
        if not args.fixedCenters:
            print '(x,y):\t(%.3f, %.3f)'%(theta[0], theta[1])
        for i in xrange(1,2):
            printParams = (i, theta[(not args.fixedCenters)+i], i, theta[(not args.fixedCenters)+i+1], i,theta[(not args.fixedCenters)+i+2], i, theta[(not args.fixedCenters)+i+3])
            print 'A%d:\t%.3f\nVx%d:\t%.3f\nVy%d:\t%.3f\nP%d:\t%.3f'%printParams
        print'\n'+'--'*20
        BEs.append(be)
        prim_fits.append(prim_fit)

        #Necessary despite centers existing because because imageObj.center is in the opposite order of numpy slicing.
        if args.fixedCenters:
            c = c_y, c_x
        else:
            c = theta[1], theta[0]

        #TODO Delete, only for testing

        from matplotlib import pyplot as plt
        imPlots = []
        fig = plt.figure(figsize = (30,20))
        minVal, maxVal = 0, 0
        plt.subplot(131)
        im = plt.imshow(imageObj[primaryBand],cmap = chosen_cmap)
        minVal = min(minVal, imageObj[primaryBand].min())
        maxVal = max(maxVal, imageObj[primaryBand].max())
        imPlots.append(im)
        plt.subplot(132)
        im = plt.imshow(prim_fit,cmap = chosen_cmap)
        minVal = min(minVal, prim_fit.min())
        maxVal = max(maxVal, prim_fit.max())
        imPlots.append(im)
        plt.subplot(133)
        calc_img = imageObj[primaryBand]-prim_fit
        im = plt.imshow(calc_img,cmap = chosen_cmap)
        minVal = min(minVal,calc_img.min())
        maxVal = max(maxVal,calc_img.max())
        imPlots.append(im)
        for image in imPlots:
            image.set_clim(minVal, maxVal)
        cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
        fig.colorbar(imPlots[0], cax, orientation='horizontal')
        fig.suptitle('%d'%1)
        plt.savefig(outputdir+imageObj.imageID+'_%d_'%1+'fullModel.png')
        plt.show()

        #estimate how much "signal" we have, so we don't overfit
        #area enclosed within nSigma
        pixelsPerParam = 10
        nSigma = 2

        varX, varY = theta[1+2*(not args.fixedCenters):3+2*(not args.fixedCenters)]
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
        #TODO delete
        #maxGaussians = 4
        for n in xrange(2,maxGaussians+1):
            prim_fit, theta, be = mcmcFit(imageObj[primaryBand], n, c_x, c_y,not args.fixedCenters, filename = name)

            print 'Gaussian #%d'%n
            if not args.fixedCenters:
                print '(x,y):\t(%.3f, %.3f)'%(theta[0], theta[1])
            for i in xrange(1,n+1):
                printParams = (i, theta[(not args.fixedCenters)+i], i, theta[(not args.fixedCenters)+i+1*n], i,theta[(not args.fixedCenters)+i+2*n], i, theta[(not args.fixedCenters)+i+3*n])
                print 'A%d:\t%.3f\nVx%d:\t%.3f\nVy%d:\t%.3f\nP%d:\t%.3f'%printParams
            print'\n'+'--'*20

            BEs.append(be)
            prim_fits.append(prim_fit)

            #TODO Delete, only for testing

            from matplotlib import pyplot as plt
            imPlots = []
            fig = plt.figure(figsize = (30,20))
            minVal, maxVal = 0, 0
            plt.subplot(131)
            im = plt.imshow(imageObj[primaryBand],cmap = chosen_cmap)
            minVal = min(minVal, imageObj[primaryBand].min())
            maxVal = max(maxVal, imageObj[primaryBand].max())
            imPlots.append(im)
            plt.subplot(132)
            im = plt.imshow(prim_fit,cmap = chosen_cmap)
            minVal = min(minVal, prim_fit.min())
            maxVal = max(maxVal, prim_fit.max())
            imPlots.append(im)
            plt.subplot(133)
            calc_img = imageObj[primaryBand]-prim_fit
            im = plt.imshow(calc_img,cmap = chosen_cmap)
            minVal = min(minVal,calc_img.min())
            maxVal = max(maxVal,calc_img.max())
            imPlots.append(im)
            for image in imPlots:
                image.set_clim(minVal, maxVal)
            cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
            fig.colorbar(imPlots[0],cax, orientation = 'horizontal')

            fig.suptitle('%d'%n)
            plt.savefig(outputdir+imageObj.imageID+'_%d_'%n+'fullModel.png')
            plt.show()


            print 'Diff: %.3f\t Old: %.3f\t New: %.3f'%(BEs[-1]-BEs[-2], BEs[-1], BEs[-2])
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
        print 'Best model N = %d'%(bestArg+1)
        if secondaryBand is not None:
            prim_fit_scaled = prim_fit*imageObj[secondaryBand][c]/imageObj[primaryBand][c]
            calc_img = imageObj[secondaryBand] - prim_fit_scaled
            calcImgDict[secondaryBand] = calc_img
        '''
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
        '''

        if args.subtraction:
            from matplotlib import pyplot as plt
            for band in bands:
                plt.figure()
                im = plt.imshow(calcImgDict[band],cmap = chosen_cmap)
                plt.colorbar(im)
                if not args.fixedCenters:
                    plt.scatter(theta[0], theta[1], color = 'm')

                #If leaving at this spot, fix filename
                plt.savefig(outputdir+imageObj.imageID+'_'+band+'_subtraction.png')
                plt.show()
                plt.clf()
                plt.close()


        if args.subtractionData:
            import numpy as np
            np.savetxt(outputdir+imageObj.imageID+'_'+band+'_residualData', calc_img)

        #TODO Plotting functionality here
        #TODO Have this flagged on/off. We don't need to check for lenses on all of em.
        #check for lens properties
        from goodnessOfFit import goodnessOfFit
        goodnessOfFit(prim_fit, imageObj[primaryBand], 4*(bestArg+1)+2*(not args.fixedCenters), 1)
        lens = residualID(calc_img, c[1], c[0])
        print 'The image %s represents a lens: %s'%(imageObj.imageID, str(lens))

if __name__ == '__main__':
    main()
