#! /usr/bin/env python
#@Author Sean McLaughlin

desc = '''
This module crops the image around the center so that only the object of interest is contained.
It contains the function cropImage, which takes a numpy array and the center coordinates as input.

cropImage(image, c_x, c_y)

Returns a numpy array of a cropped image, and the readjusted c_x and c_y

There is an optional plot functionality that saves the cropped image to file.

cropImage(image, c_x, c_y, plot = True, filename = fname)

This module can be run as main. In this case it takes a filename for a fits file and the center coordinates, as well
as a filename where to save the plot.

python cropImage.py filename c_x c_y savename

An optional parameter in both cases in the sideLength, which is the size of the cropped image. 

cropImage(image, c_x, c_y, sideLength = 20)

python cropimage.py filename c_x c_y savename 20

If none is given, it will be set to 30.
'''

import numpy as np
from matplotlib import pyplot as plt

def cropImage(image, c_x, c_y, sideLength= 30, plot = False, filename = None):

    img_y, img_x = image.shape

    xLims = [int(c_x - .5*sideLength),int( c_x + .5*sideLength)]  
    yLims = [int(c_y - .5*sideLength), int(c_y + .5*sideLength)]  

    #if any of the lims are out of bounds, make them the edge instead
    for i in xrange(2):
        if xLims[i] <0:
            xLims[i] = 0
        if xLims[i]>img_x:
            xLims[i] = img_x
        if yLims[i] <0:
            yLims[i] = 0
        if yLims[i]>img_y:
            yLims[i] = img_y
                                                                    
    image = image[yLims[0]:yLims[1], xLims[0]:xLims[1]]                
    #readujust the centers. They've moved now that the image has been cropped
    c_y, c_x = c_y-yLims[0], c_x - xLims[0]

    if plot:
        plt.figure()
        im = plt.imshow(image)
        plt.colorbar(im)
        plt.scatter(c_x, c_y, color = 'k')
        plt.savefig(filename)
        plt.clf()
        plt.close()

    return image, c_x, c_y

if __name__ == '__main__':
    import pyfits
    import argparse

    #process inputs cleanly
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('filename', metavar = 'fname', type = str, help = 'The name of the fits file')
    parser.add_argument('center_x', metavar = 'c_x', type = int, help = 'The center in x')
    parser.add_argument('center_y', metavar = 'c_y', type = int, help = 'The center in y')
    parser.add_argument('savename', metavar = 'sname', type = str, help = 'The name of the file to save the plot to.')
    parser.add_argument('sideLength', metavar = 'sLength', nargs = '?',type = int, help = 'The optional length of a side for the crop', default = 30)

    args = parser.parse_args()

    filename = args.filename

    try:
        fitsImage = pyfits.open(filename)
    except IOError:
        print 'ERROR: Invalid filename.'
        from sys import exit
        exit(-1) 

    #TODO add default for cetners, and if none given call findCenter.py

    image = fitsImage[0].data
    c_y, c_x = args.center_y, args.center_x

    image, c_x, c_y = cropImage(image, c_x, c_y, args.sideLength, True, args.savename) 
