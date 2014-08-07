#makes a fits image that is a series of Gaussians on a black field.
#I'm going to have one "main" galaxy with another Gaussian a small bit away.  

import pyfits
import os
import numpy as np

if 'lens.fits' in os.listdir(os.getcwd()):
	os.remove('lens.fits')

def gaussian(x,y):
	return A1*np.exp(-1./(R1**2)*((x-CENTER[0])**2+(y-CENTER[1])**2)) \
			+ A2*np.exp(-1./(R2**2)*((x-CENTER2[0])**2+(y-CENTER2[1])**2)) 

IMAGE_SIZE = 200
CENTER = (IMAGE_SIZE/2, IMAGE_SIZE/2)
CENTER2 = (CENTER[0]+IMAGE_SIZE/10., CENTER[1]+IMAGE_SIZE/10.)
A1 = 10
R1 = 25 
A2 = 3
R2 = 10 

image = [ [ gaussian(i,j) for i in xrange(IMAGE_SIZE) ] for j in xrange(IMAGE_SIZE)]

image = np.array(image)

noise = np.random.normal(size = image.shape)

image+=noise

hdu = pyfits.PrimaryHDU(image)
hdu.writeto('lens.fits')
