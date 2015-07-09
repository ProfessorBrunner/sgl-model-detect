#! /usr/bin/env python
#@Author Sean McLaughlin

gooeyON = False #Wheter or not to use the Gooey module or stick with standard CLI
#In addition to setting as True, make sure to uncomment the decorator in line 11!

if gooeyON:
    from gooey import Gooey, GooeyParser

#Gooey doesn't seem to let the columns change...
#@Gooey( required_cols= 3)
def main():
    desc = '''
    This program models the light profile of a galaxy in a fits image using a Mixture-of-Gaussian model. It can use either MCMC or non-linear least squares.
    '''
    #Uncomment for a GUI-less version!
    if gooeyON:
        parser = GooeyParser(description = desc)
    else:
        from argparse import ArgumentParser
        parser = ArgumentParser(description = desc)

    if gooeyON:#very similar to the segment below, only difference is in the GUI it adds file choosers.
        parser.add_argument('filename', metavar = 'filename', widget = 'FileChooser', type = str,help = \
                                        'Either a fits filename or a directory of files.\n NOTE: The GUI will not allow you to choose a directory. Choose a file and edit the text.')

        parser.add_argument('centers', metavar = 'centers', widget = 'FileChooser',type = str,help = \
                                        'Either a filename or a comma separate pair of coordinates for x,y.')

        parser.add_argument('output', metavar = 'output', widget = 'DirChooser',  type = str, help = \
                                        'Location to store the programs output.')

    else:

        parser.add_argument('filename', metavar = 'filename',  type = str,help = \
                                        'Either a fits filename or a directory of files.\n NOTE: The GUI will not allow you to choose a directory. Choose a file and edit the text.')

        parser.add_argument('centers', metavar = 'centers', type = str,help = \
                                        'Either a filename or a comma separate pair of coordinates for x,y.')

        parser.add_argument('output', metavar = 'output',  type = str, help = \
                                        'Location to store the programs output.')

    parser.add_argument('modeler', metavar = 'modeler', type = str, choices = ['MCMC', 'NLSQ'], help = \
                                    'Modeler used to perform fit. Either MCMC or NLSQ.')

    walkersChoices = xrange(0,4500,500)
    #GUI requires strings, CLI requires ints
    if gooeyON:
        walkersChoices = [ str(i) for i in walkersChoices]
    parser.add_argument('n_walkers', metavar = 'n_walkers', choices = walkersChoices, type = int, help =\
                                    'Number of walkers to use in the MCMC sampler. If NLSQ is chosen choice does not matter. Multiples of 500 up to 4000 allowed.')

    parser.add_argument('format', metavar = 'format', type = str, choices = ['CFHTLS', 'SDSS', 'Toy'], help = \
                                    'The format of images to use. Choice of CFHTLS, SDSS or Toy')

    parser.add_argument('primaryBand', metavar = 'primaryBand',  choices = ['u','g','r','i', 'z', 't'], type = str, help = \
                                    "What band to use in the fit. Choice of u,g,r,i,z, or t")
    #NOTE Not sure how the GUI handles the optional args
    parser.add_argument('secondaryBand', metavar = 'secondaryBand',  choices = ['None','u','g','r','i', 'z', 't'], type = str, nargs = '?', default = 'None', help = \
                                    "Optional. If option other than 'None' chosen, will subtract a scaled version of the model in primary band from this one and search for lens residuls.")


    parser.add_argument('--fixedCenters', dest = 'fixedCenters', action = 'store_true',help =\
                                    'Use the given centers as fixed values. Otherwise, they will be treated as initial guesses')
    parser.add_argument('--showImages', dest = 'showImages', action = 'store_true', help = \
                                    'Show images rather than save them to file.')
    parser.add_argument('--cutout', dest = 'cutout', action = 'store_true',\
                         help = 'Save a .png of the original cutout to file.')
    parser.add_argument('--cutoutData', dest = 'cutoutData', action = 'store_true',\
                         help = 'Save the raw data of the cutout to file.')
    parser.add_argument('--subtraction', dest = 'subtraction', action = 'store_true',\
                         help = 'Store a .png of the model subtraction from the original image.')
    parser.add_argument('--subtractionData', dest = 'subtractionData', action = 'store_true',\
                                help = 'Save the raw image data to file of the residuals.')
    parser.add_argument('--chain', dest = 'chain', action = 'store_true',\
                     help = 'Store the markov chain in the output directory. OCCUPIES A LOT OF FILESPACE.')

    args = parser.parse_args()
    SHOW_IMAGES = args.showImages

    if not SHOW_IMAGES:
        import matplotlib as mpl
        mpl.use('Agg')

    from mcmcFit import mcmcFit, printTheta
    from nlsqFit import nlsqFit, goodnessOfFit
    from residualID import residualID
    import imageClass
    import numpy as np
    from matplotlib import pyplot as plt
    import seaborn
    seaborn.set()
    import os
    from time import time

    chosen_cmap = 'gnuplot2'

    def plotSingleImage(imageDict, bands, chosen_cmap, outputdir,id, name, show = False):
        'Plot one image. Used for both subtraction and cutout.'
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

    #TODO Make the font larger and add titles to the subplots
    def plotFullModel(rawImage, model, N, id, outputdir, chosen_cmap, show = False ):
        'Function that creates a more elaborate plot of the image, model, and residuals.'
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

        #Same min and max for all images
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
            maxGaussians = 2 #minimum value of 2
        elif maxGaussians > 10:
            maxGaussians = 10 #upper value of 10. Arbitrarily chosen and can be extended.
        print 'Max Gaussians = %d'%maxGaussians

        return maxGaussians

    filename = args.filename
    outputdir = args.output

    files = [filename, outputdir]

    primaryBand = args.primaryBand
    secondaryBand = None if args.secondaryBand == 'None' else args.secondaryBand
    bands = [primaryBand]
    if secondaryBand is not None:
        bands.append(secondaryBand)

    #determine which center finding method to use
    isCoordinates = args.centers.find(',') != -1

    if isCoordinates:
        splitCenters = args.centers.split(',')
        c_x = int(splitCenters[0].strip())
        c_y = int(splitCenters[1].strip())

    #Then we've been given a list of centers that will have to be read in.
    else:
        files.append(args.centers)
        centers = args.centers

    #Directories need a slash at the end of them
    if os.path.isdir(filename) and filename[-1] != '/':
        filename+='/'

    #check files for existance
    for f in files:
        if not os.path.exists(f):
            print 'ERROR: Invalied path %s'%f
            from sys import exit
            exit(-1)

    #Make a dictionary of the Galaxy's image coordinates.
    if not isCoordinates :
        centerDict = {}
        idx = 2 if args.format == 'CFHTLS' else 1 #The catalogs are formatted slightly differently
        # This was easier than fixing the catalogs, but a standard format may be preferable.
        with open(centers) as f:
            for line in f:
                splitLine = line.strip().split(' ')
                coords = [float(x) for x in splitLine[idx:]]
                centerDict[splitLine[0]] = coords
    else:
        centerDict = None

    #both functions have the same API, so they can be swapped out no problem.
    modelerDict = {'MCMC': mcmcFit, 'NLSQ': nlsqFit}
    fitter = modelerDict[args.modeler]

    imageClassDict = {'CFHTLS': imageClass.CFHTLS, 'SDSS': imageClass.SDSS, 'Toy': imageClass.Toy}
    #the appropriate formatting for these objects
    imgClass = imageClassDict[args.format]
    imageDict = {}

    #load in filenames from directory
    if os.path.isdir(filename):
        #if the input is a directory
        fileList = os.listdir(filename)
        fileDirectory = filename
        for fname in fileList:
            #Only look at .fits images
            if fname[-5:] != '.fits':
                continue
            obj = imgClass(fileDirectory+fname)

            if (primaryBand not in obj.images and secondaryBand not in obj.images):
                #it's got neither of the bands we're interested in, so skip it!
                continue

            if obj.imageID in imageDict:
                #We already loaded a different band of this object
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
        print '-'*30
        print '*'*30
        print '_'*30
        print 'Image ID : %s'%imageObj.imageID

        #make a folder for all outputs
        #All outputs for this object will go in the same folder!
        imOutputDir = outputdir+imageObj.imageID+'/'
        if os.path.isdir(imOutputDir):
            for f in os.listdir(imOutputDir):
                os.remove(imOutputDir+f) #Delete all the stuff that's in there first
                #Maybe have some sort of warning first?
        else:
            os.mkdir(imOutputDir)

        if isCoordinates:
            coords = (c_x, c_y)
        else:
            coords = None

        #adjusts the center estimate, or finds one outright if none were given.
        try:
            imageObj.calculateCenter(coords,centerDict)
        except KeyError:
            print 'No center for %s'%imageObj.imageID
            continue

        #crop the image down to size
        imageObj.cropImage()

        #Plot the Requested Cutout
        if args.cutout:
            plotSingleImage(imageObj, bands, chosen_cmap, imOutputDir,imageObj.imageID, 'cutout', show = SHOW_IMAGES)

        if args.cutoutData:
            for band in bands:
                np.savetxt(imOutputDir+imageObj.imageID+'_'+band+'_cutoutData', imageObj.images[band])

        stats = []
        prim_fits = []
        thetas = []
        #perform first fit with 1 Gaussian
        #Then, use to charecterize max number of parameters
        c_x, c_y = imageObj.center #Note that this will have changed after the crop
        print 'Fitting now'
        prim_fit,theta, stat  = fitter(imageObj[primaryBand], 1, c_x, c_y, not args.fixedCenters,n_walkers = args.n_walkers, dirname = imOutputDir, id = imageObj.imageID, chain = args.chain)

        printTheta(theta, 1, movingCenters= not args.fixedCenters)

        stats.append(stat)
        prim_fits.append(prim_fit)
        thetas.append(theta)

        plotFullModel(imageObj[primaryBand], prim_fit, 1, imageObj.imageID, imOutputDir, chosen_cmap, show = SHOW_IMAGES)

        maxGaussians = caculateMaxGaussians(theta)

        for n in xrange(2,maxGaussians+1):
            prim_fit, theta, stat = fitter(imageObj[primaryBand], n, c_x, c_y,not args.fixedCenters,n_walkers = args.n_walkers, dirname = imOutputDir, id = imageObj.imageID, chain = args.chain)

            printTheta(theta,2, movingCenters= not args.fixedCenters)

            if np.isnan(stat):
                print 'NaN raised for Gaussian %d'%n
                break

            stats.append(stat)
            prim_fits.append(prim_fit)
            thetas.append(theta)

            plotFullModel(imageObj[primaryBand], prim_fit, n, imageObj.imageID, imOutputDir, chosen_cmap, show = SHOW_IMAGES)

            print 'Old BE: %.3f\t New BE: %.3f'%(stats[-2], stats[-1])

            if args.modeler == 'NLSQ':
                #do a chi2 test
                test = abs(stats[-1]-1)> abs(stats[-2]-1)
            else: #Bayes Factor
                test = stats[-1]<stats[-2] #equivalent to BF < 1

            if test: #new Model is worse!
                break
                #pass

        if args.modeler == 'NLSQ':
            bestArg = np.argmin(np.abs(np.array(stats)-1))
        else:
            bestArg = np.argmax(np.array(stats))

        bestTheta = thetas[bestArg]
        #TODO Instead of writing a numpy array, perhaps using printTheta to do a more visual print would be better.
        np.savetxt(imOutputDir+imageObj.imageID+'_'+primaryBand+'_theta', bestTheta, delimiter=',')

        prim_fit = prim_fits[bestArg]
        print 'Best model N = %d'%(bestArg+1)

        calcImgDict = {}
        calc_img = imageObj[primaryBand] - prim_fit
        calcImgDict[primaryBand] = calc_img

        if secondaryBand is not None:

            if args.fixedCenters:
                c = (imageObj.center[1], imageObj.center[0])
            else:
                c = int(bestTheta[1]), int(bestTheta[0])

            prim_fit_scaled = prim_fit*imageObj[secondaryBand][c]/imageObj[primaryBand][c]
            calc_img = imageObj[secondaryBand] - prim_fit_scaled
            calcImgDict[secondaryBand] = calc_img

        if args.subtraction:
            plotSingleImage(calcImgDict, bands, chosen_cmap, imOutputDir, 'subtraction',imageObj.imageID, show = SHOW_IMAGES)

        if args.subtractionData:
            np.savetxt(imOutputDir+imageObj.imageID+'_'+primaryBand+'_residualData', calc_img)

        goodnessOfFit(prim_fit, imageObj[primaryBand], 4*(bestArg+1)+2*(not args.fixedCenters), 1)

        #TODO Plotting functionality here
        #check for lens properties
        if secondaryBand is not None:
            lens = residualID(calc_img, c[1], c[0])
            print 'The image %s represents a lens: %s'%(imageObj.imageID, str(lens))

if __name__ == '__main__':
    main()
