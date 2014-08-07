#makes a fits image that is a white circle in a black field

import pyfits
import os
import numpy as np

if 'circle.fits' in os.listdir(os.getcwd()):
	os.remove('circle.fits')

CIRCLE_RAD = 8 
IMAGE_SIZE = 200
CENTER = (IMAGE_SIZE/2, IMAGE_SIZE/2)

def in_circle(x,y): #determines if x,y is in the central circle
	return int((x-CENTER[0])**2+(y-CENTER[1])**2 <=IMAGE_SIZE*CIRCLE_RAD)

image = [ [ in_circle(i,j) for i in xrange(IMAGE_SIZE) ] for j in xrange(IMAGE_SIZE)]

image = np.array(image)

hdu = pyfits.PrimaryHDU(image)
hdu.writeto('circle.fits')

#now i'm adding functionality to make a "sphere" of brightness
if 'sphere.fits' in os.listdir(os.getcwd()):
	os.remove('sphere.fits')

def in_sphere(x,y): 
	val =-1*int((x-CENTER[0])**2+(y-CENTER[1])**2-IMAGE_SIZE*CIRCLE_RAD)
	return val if val>0 else 0

image = [ [in_sphere(i,j) for i in xrange(IMAGE_SIZE)] for j in xrange(IMAGE_SIZE)] 
image = np.array(image)
hdu = pyfits.PrimaryHDU(image)
hdu.writeto('sphere.fits')
