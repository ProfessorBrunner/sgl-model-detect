#! /usr/bin/env python
#@Author Sean McLaughlin
#This file contains the Image class, which contains all the details of loading in fits images
#It will also include a few subclasses depending on filename syntax

import pyfits
from findCenter import findCenter
from cropImage import cropImage

class Image(object):
    
    def __init__(self, filename):

        self.filenames = [filename]

        self.imageID = self._getImageID(filename) 
        band = self._getBand(filename) 
        fitsImage = pyfits.open(filename)
        self.images = { band : fitsImage[0].data}

    def __getitem__(self, key):
        return self.images[key]

    def __setitem__(self, key, value):
        self.images[key] = value

    def __delitem__(self, key):
        del self.images[key]

    def addImage(self,filename):
        self.filenames.append(filename)

        band = self._getBand(filename) 
        fitsImage = pyfits.open(filename)
        self.images[band] = fitsImage[0].data

    def calculateCenters(self, coords = None, galaxyDict = None):
        
        image = self.images.values()[0] #get first image
        #Is there a way I should make use of the multiple images?
        if coords is not None:
            c_x, c_y = coords
        elif galaxyDict is not None:
            c_x, c_y = inputDict['galaxyDict'][self.imageID]
        else:
            c_y, c_x = findCenter(image)

        self.center = (c_x, c_y)

    def cropImage(self, plot = False, output = ''):
        for band, image in self.images.iteritems():
            image, c_x, c_y = cropImage(image, self.center[0], self.center[1], plot = plot, filename = output + self.imageID+'_'+band+'_cutout.png')
            self.images[band] = image
        self.center = (c_x, c_y)

    def _getImageID(self, filename):
        #Implemented in subclasses
        return None
        
    def _getBand(self, filename):
        #Implemented in subclasses
        return '' 

class CFHTLS(Image):
   
    def _getImageID(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex], filename[lineIndex:]
        return baseName[:-7]
    
    def _getBand(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex], filename[lineIndex:]
        return baseName[-6]
        
class SDSS(Image):

    def _getImageID(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex], filename[lineIndex:]
        return baseName[9:-5]
        
    def _getBand(self, filename):
        lineIndex = filename.rfind('/')
        fileDirectory, baseName= filename[:lineIndex], filename[lineIndex:]
        return baseName[7]
