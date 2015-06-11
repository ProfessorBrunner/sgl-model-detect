#! /usr/bin/env python
#@Author Sean McLaughlin
#This file contains the Image class, which contains all the details of loading in fits images
#It will also include a few subclasses depending on filename syntax

import pyfits
from cropImage import cropImage

class Image(object):
    
    def __init__(self, filename):

        self.imageID = self._getImageID(filename) 
        band = self._getBand(filename)
        self.filenames = {band : filename}
        fitsImage = pyfits.open(filename)
        self.images = { band : fitsImage[0].data}

    def __getitem__(self, key):
        return self.images[key]

    def __setitem__(self, key, value):
        self.images[key] = value

    def __delitem__(self, key):
        del self.images[key]

    def addImage(self,filename):

        band = self._getBand(filename)
        self.filenames[band] = filename
        fitsImage = pyfits.open(filename)
        self.images[band] = fitsImage[0].data

    def calculateCenters(self, coords = None, galaxyDict = None):
        import numpy as np
        image = self.images.values()[0] #get first image
        #Is there a way I should make use of the multiple images?
        if coords is not None:
            c_x, c_y = coords
        elif galaxyDict is not None:
            c_x, c_y = galaxyDict[self.imageID]
        else:
            c_y, c_x = np.unravel_index(image.argmax(), image.shape) #center is just the brightest spot
            #c_y, c_x = findCenter(image)

        #sometimes the center is not exactly accurate. This part finds the maximum in the region around the center.
        dy, dx = np.unravel_index(image[c_y-2:c_y+3, c_x-2:c_x+3].argmax(), (5,5))
        dy,dx = dy-2, dx-2 #so if c_x c_y is the max dx, dy will be 0,0
        c_y, c_x = c_y+dy, c_x+dx

        self.center = (c_x, c_y)

    def cropImage(self, plot = False, output = ''):
        for band, image in self.images.iteritems():
            image, c_x, c_y = cropImage(image, self.center[0], self.center[1], plot = plot, filename = output + self.imageID+'_'+band+'_cutout.png')
            self.images[band] = image
        self.center = (c_x, c_y)

    def getOtherBand(self,bands):
        'Attempts to get other images of the same ID in different bands. Assumes only one image is loaded into the object so far'
        #Implemented in subclasses
        pass

    def _getImageID(self, filename):
        #Implemented in subclasses
        return None
        
    def _getBand(self, filename):
        #Implemented in subclasses
        return ''

class CFHTLS(Image):

    def getOtherBand(self, bands):
        'Attempts to get other images of the same ID in different bands. Assumes only one image is loaded into the object so far'
        #Note that band is a string of all bands used in the image. The program first checks which one we already have!

        for band in bands:
            if band not in self.filenames:
                break

        filename = self.filenames.values()[0]
        lineIndex = filename.rfind('/')
        fileDirectory, baseName = filename[:lineIndex+1], filename[lineIndex+1:]
        otherFilename = ''.join([fileDirectory,'CFHTLS_',self.imageID,'_', band,'.fits'])
        self.addImage(otherFilename)
   
    def _getImageID(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex+1], filename[lineIndex+1:]
        return baseName[7:-7]
    
    def _getBand(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex], filename[lineIndex:]
        return baseName[-6]
        
class SDSS(Image):

    def getOtherBand(self,bands):
        'Attempts to get other images of the same ID in different bands. Assumes only one image is loaded into the object so far'
        #Note that band is a string of all bands used in the image. The program first checks which one we already have!

        for band in bands:
            if band not in self.filenames:
                break

        filename = self.filenames.values()[0]
        lineIndex = filename.rfind('/')
        fileDirectory, baseName = filename[:lineIndex+1], filename[lineIndex+1:]
        otherFilename = ''.join([fileDirectory,'frame-',band,'-',self.imageID,'.fits'])
        self.addImage(otherFilename)

    def _getImageID(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex], filename[lineIndex:]
        return baseName[9:-5]
        
    def _getBand(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex], filename[lineIndex:]
        return baseName[7]
