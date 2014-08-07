#! /usr/bin/env python
#@Author Sean McLaughlin
desc ='''
The standard way to find the center of an astronomical image is with moments.
However, this doesn't work very well if there are multiple objects in the image.
This code finds the center of the brightest object in the image by finding the moments of 
smaller and smmaller boxes centered around the previous moment.

The funcion findCenter can be imported.

findCenter(image) 

Where image is a numpy array. It contains an optional plot arguement, where the plot is saved to filename.

findCenter(image, plot = True, filename = 'center.png')

This module can also be run as main. It automatically plots in this case.

python findCenter.py  imageFilename saveFilename
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm as lm
from matplotlib import cm

def findCenter(image, plot = False, filename = None):

        if plot and filename is None:
                raise ValueError('Plot is true but not filename given in findCenter.')

        img_y, img_x = image.shape
	#start at the maximum position
	maxPos =  np.unravel_index(image.argmax(), image.shape)
	boxRatio = .5 #the box around the max will be half the size of the image.

	xLims = (maxPos[1]-img_x*boxRatio*.5, maxPos[1]+img_x*boxRatio*.5)
	#rescale so the limits do not expand beyond image boundaries.
	xLims = (int(xLims[0]) if xLims[0]>=0 else 0, int(xLims[1]) if xLims[1]<=img_x else img_x)

	yLims = (maxPos[0]-img_y*boxRatio*.5, maxPos[0]+img_y*boxRatio*.5)
	yLims = (int(yLims[0]) if yLims[0]>=0 else 0, int(yLims[1]) if yLims[1]<=img_y else img_y)

	threshold = image.mean()+image.std()
	N = 10 #10 iterations
	
	if plot:
		fig = plt.figure()
		plt.scatter(maxPos[1], maxPos[0], c= 'b')

	for n in xrange(N):
		xmean = ymean = totval = 0

		for i in xrange(*xLims):
		    for j in xrange(*yLims):
		        if image[j,i]>threshold:
                            val = image[j,i]
		            xmean+=i*val
		            ymean+=j*val
		            totval +=val


		xmean, ymean = xmean/float(totval), ymean/float(totval)

		if plot:
			plt.imshow(image,\
				    norm = lm(image.mean() + .5*image.std(), image.max(), clip = 'True'),\
				    cmap = cm.gray, origin = 'lower')
			plt.scatter(xmean, ymean, c='r')
			plt.plot([xLims[0] for i in xrange(100)], np.linspace(yLims[0],yLims[1], 100), color = 'm', linewidth= 2, alpha = .2)
			plt.plot([xLims[1] for i in xrange(100)], np.linspace(yLims[0], yLims[1], 100), color ='m', linewidth= 2, alpha = .2)
			plt.plot(np.linspace(xLims[0],xLims[1], 100),[yLims[0] for i in xrange(100)], color = 'm', linewidth= 2, alpha = .2)
			plt.plot(np.linspace(xLims[0],xLims[1], 100),[yLims[1] for i in xrange(100)], color = 'm', linewidth= 2, alpha = .2)

		maxPos = (ymean, xmean) 
		boxRatio = .8*boxRatio #the box around the max will be half the size of the image.

		xLims = (maxPos[1]-img_x*boxRatio*.5, maxPos[1]+img_x*boxRatio*.5)
		xLims = (int(xLims[0]) if xLims[0]>=0 else 0, int(xLims[1]) if xLims[1]<=img_x else img_x)

		yLims = (maxPos[0]-img_y*boxRatio*.5, maxPos[0]+img_y*boxRatio*.5)
		yLims = (int(yLims[0]) if yLims[0]>=0 else 0, int(yLims[1]) if yLims[1]<=img_y else img_y)
	
	if plot:
		#plt.show()
		plt.savefig(filename)
		
	return maxPos

if __name__ == '__main__':
	import pyfits
	import argparse
    #process inputs cleanly
	parser = argparse.ArgumentParser(description = desc)
	parser.add_argument('filename', metavar = 'fname', type = str, help = 'The name of the fits file')
    parser.add_argument('savename', metavar = 'sname', type = str, help = 'The name of the file to save the plot to.')

	args = parser.parse_args()
	fname = args.filename
	if fname.split('.')[-1] != 'fits':
	#if they are not fits files
		from sys import exit
		print "ERROR: %s not a valid fits filename."%fname 
		exit(-1)

	try:
		fitsImage = pyfits.open(fname)
	except IOError:
		print 'ERROR: Invalid filename: %s'%fname
		from sys import exit
		exit(-1)

	image = fitsImage[0].data
	find_Center(image, plot= True, args.savename)
        
