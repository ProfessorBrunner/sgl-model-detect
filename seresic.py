import numpy as np
import scipy as sp
import scipy.ndimage
from scipy import misc
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import sys
from dcor import * 

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Experimental algorithm for finding the center of a foreground
lensing galaxy.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''

def plot_polar_image(data, origin=None):
	#Plots an image reprojected into polar coordinates

	polar_grid, r, theta = reproject_image_into_polar(data, origin)
	plt.figure()
	plt.imshow(polar_grid,extent=(theta.min(), theta.max(), r.max(), r.min()))
	plt.axis('auto')
	plt.xlabel('Theta Coordinate (radians)')
	plt.ylabel('R Coordinate (pixels)')
	plt.title('Image in Polar Coordinates')
	return polar_grid

def index_coords(data, origin=None):
	#Creates x & y coords for the indicies in a numpy array
	ny, nx = data.shape[:2]
	if origin is None:
		origin_x, origin_y = nx // 2, ny // 2
	else:
		origin_x, origin_y = origin
	x, y = np.meshgrid(np.arange(nx), np.arange(ny))
	x -= origin_x
	y -= origin_y
	return x, y

def cart2polar(x, y):
	r = np.sqrt(x**2 + y**2)
	theta = np.arctan2(y, x)
	return r, theta

def polar2cart(r, theta):
	x = r * np.cos(theta)
	y = r * np.sin(theta)
	return x, y


def reproject_image_into_polar(data, origin=None):
	#Reprojects a 3D numpy array ("data") into a polar coordinate system.
	ny, nx = data.shape[:2]

	if origin is None:
		origin = (nx//2, ny//2)
    																	# Determine range of r and theta
	x, y = index_coords(data, origin=origin)
	r, theta = cart2polar(x, y)

	r_i = np.linspace(r.min(), r.max(), nx)								# Making a grid
	theta_i = np.linspace(theta.min(), theta.max(), ny)
	theta_grid, r_grid = np.meshgrid(theta_i, r_i)

    																	# Project the r and theta grid back into pixel coordinates
	xi, yi = polar2cart(r_grid, theta_grid)
	xi += origin[0] 													# Shift the origin back 
	yi += origin[1] 
	xi, yi = xi.flatten(), yi.flatten()
	coords = np.vstack((xi, yi))

	zi = sp.ndimage.map_coordinates(data, coords, order=1).reshape((nx, ny))
	return zi, r_i, theta_i

def sersic(x,k,n,s):

	return s * math.e ** (-1.0 * k * (  x ** (1.0/n) ) ) 

def get_statistics(data):
	m , n = data.shape
	row_mean = []
	row_median = []
	row_index = []
	for i in range (0,m-1):
		row_mean.append(np.mean(data[i,:]))
		row_median.append(np.median(data[i,:]))
		row_index.append(i)
	return row_mean, row_median, row_index

def fit_profile(row_mean,row_index):
	popt, pcov = curve_fit(sersic, np.asarray(row_index), np.asarray(row_mean))
	Model_Fit =[]
	for i in range (0,len(row_index)):
		Model_Fit.append(sersic(row_index[i],popt[0],popt[1],popt[2])) 
	return Model_Fit
	

if __name__ == '__main__':

	if (len(sys.argv)!=2):
		print 'Please specify one input file'
		exit(0)
	
	z = misc.imread(sys.argv[1],True)
	m , n = z.shape

	Corr = []
	cVal = 0

	Max = 0
	Max_x = 0
	Max_y = 0 
	Val_sum = 0
	for i in range (364,365):
		for j in range (413,n-1):
			print 'Coordinates: ' + str( i ) + ', ' +str (j)
			Val_sum += z[i][j]												#Find Sum to calculate Mean
			try:			
				polar_image = plot_polar_image(z, (i,j))
				row_mean, row_median, row_index = get_statistics(polar_image)
			
				Model_Fit = fit_profile(row_mean,row_index)
				cVal = dcor(np.asarray(Model_Fit),np.asarray(row_mean)).find_correlation()
				Corr.append([cVal,i,j])				
				print cVal
				if(cVal>Max):
					Max = cVal
					Max_x = i
					Max_y = j 
			except(RuntimeError):
				print 'Optimized parameters not found'				
				pass
				

	Mean =  Val_sum / (float)( m * n )
	print Corr
	print 'MAX'
	print m, n	
	print (Max_x, Max_y, cVal)	
	polar_image = plot_polar_image(z, (Max_x,Max_y))
	row_mean, row_median, row_index = get_statistics(polar_image)
	Model_Fit = fit_profile(row_mean,row_index)



	'''
	Plotting
	'''
	Residual = []
	for i in range (0,len(row_mean)):
		Residual.append(abs(Model_Fit[i] - row_mean[i]))

	plt.figure()
	plt.plot(row_index,row_mean,label='Mean Row Intensity')
	plt.plot(row_index,row_median,label='Median Row Intensity')
	plt.plot(row_index,Model_Fit,label='Sersic Profile')
	plt.plot(row_index,Residual,label='Residual Profile')
	plt.legend()
	plt.xlabel('Row Number')
	plt.ylabel('Intensity')
	plt.show()
