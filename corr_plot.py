import sys
import pylab
from scipy import interpolate
import random
from dcor import * 
import operator

'''
Class for light curve segment
'''
class Segment(object):
	def __init__(self, intensity_function, intensity_data, time_data):
		self.intensity_function = intensity_function
		self.time_start  = min(time_data)
		self.time_end 	 = max(time_data)
		self.time_data   = time_data
		self.intensity	 = intensity_data

	'''
	Shifts a light curve segment
	'''	
	def shift_segment(self, days):
		self.time_start  = self.time_start + days
		self.time_end 	 = self.time_end + days
		
		for x in range(0,len(self.time_data)):
			self.time_data[x] = self.time_data[x] + days

		self.intensity_function = interpolate.UnivariateSpline(self.time_data, self.intensity, k = 3, s = 0.04)

	'''
	Returns the intensities corresponding to temporal overlap of light curves
	'''
	def intersection(self, segment):
		#if no overlap return empty
		if ( (self.time_end <= segment.getStartTime()) or (self.time_start >= segment.getEndTime()) ) :	
			return [] , []
		else: 
			start, end = -1, -1
			#find the common region			
			if(segment.getEndTime()<= self.time_end) and (self.time_start >= segment.getStartTime()):
				start	= self.time_start
				end		= segment.getEndTime()
			elif(segment.getEndTime()<= self.time_end) and (self.time_start <= segment.getStartTime()):
				start	= segment.getStartTime()
				end		= segment.getEndTime()
			elif(segment.getEndTime()>= self.time_end) and (self.time_start <= segment.getStartTime()):
				start	= segment.getStartTime()
				end		= self.time_end
			elif(segment.getEndTime()>= self.time_end) and (self.time_start >= segment.getStartTime()):
				start	= self.time_start
				end		= self.time_end
			
			common_segment1 = []
			common_segment2 = []	
			for t in range (int(start), int(end)):	
				common_segment1.append(self.intensity_function(t))
				common_segment2.append(segment.getIntensity(t))
			return common_segment1 , common_segment2

	'''
	Returns the intensitiy at a given time
	'''
	def getIntensity(self, time):
		return self.intensity_function(time)
	
	def getStartTime(self):
		return self.time_start

	def getEndTime(self):
		return self.time_end

	def getSegmentLength(self):
		return self.time_end - self.time_start

	def display(self):
		print 'Segment start:'
		print self.time_start
		print 'Segment end:'
		print self.time_end

'''
Returns a lists of elements having the format
[ Time, Lightcurve A, Error in A, Lightcurve B, Error in B ]
Header of input file is ignored (first two lines) 
'''
def load_data():
	if (len(sys.argv)!=2):
		print 'Please specify one input file'
		exit(0)
	time = []
	lc_a = []
	er_a = []
	lc_b = []
	er_b = []
	infile = open(sys.argv[1], 'r')
	lines = infile.readlines()
	Data = []	
	for i in range (2,len(lines)):
		elements = lines[i].split()
		time.append(float(elements[0]))
		lc_a.append(float(elements[1]))
		er_a.append(float(elements[2]))
		lc_b.append(float(elements[3]))
		er_b.append(float(elements[4]))
	infile.close()
	return time, lc_a, er_a, lc_b, er_b

'''
Plots the data points 
'''
def plot_data_points(time, lc_a, lc_b, avg_a, avg_b):
	pylab.plot(time, [l + avg_a for l in lc_a], 'rx' , label='Curve A')
	pylab.plot(time, [l + avg_b for l in lc_b], 'bx', label='Curve B')


'''
Plots a Light Curve Segment 
'''
def plot_segment(LC, avg):
	t_range = []
	intensity = []	
	for t in range(int(LC.getStartTime()),int(LC.getEndTime())): 	
		t_range.append(t)
		intensity.append(avg + LC.getIntensity(t))		
	pylab.plot(t_range, intensity,'g')

def show_plots():
	pylab.legend(loc=4)
	pylab.xlabel('Modified Heliocentric Julian Date (MHJD)')
	pylab.ylabel('Magnitude (relative)')
	pylab.title(sys.argv[1])
	pylab.show()


def find_shift(LC_1, LC_2, Max_Shift, Shift_Step):
	
	Corr_Data = []
	i = 0

	#Shift Light Curve 1		
	while (i < Max_Shift):
		print 'Shift '+ str(i) + '...'		
		
		Intersection_vector_1 = []
		Intersection_vector_2 = []
		
		#Shift curve 1 forward		
		for s in LC_1:
			s.shift_segment(Shift_Step)

		#Find intersections
		for s1 in LC_1:		
			for s2 in LC_2:
				c1, c2 = s1.intersection(s2)
				if(c1 != [] and c2 != []):
					for x in range (0,len(c1)):
						Intersection_vector_1.append(c1[x])
						Intersection_vector_2.append(c2[x])	
		print 'Overlap Vector Length = ' + str(len(Intersection_vector_1))
		Corr_Data.append([i,dcor(Intersection_vector_1,Intersection_vector_2).find_correlation()])
		
		i = i + Shift_Step 


	for s in LC_1:
		s.shift_segment(-1 * i)

	
	#Shift Light Curve 2
	i = 0
	while (i < Max_Shift):
		print 'Shift '+ str(i) + '...'		
	
		#Shift curve 2 forward		
		Intersection_vector_1 = []
		Intersection_vector_2 = []
		
		for s in LC_2:
			s.shift_segment(Shift_Step)

		#Find intersections
		for s1 in LC_1:		
			for s2 in LC_2:
				c1, c2 = s1.intersection(s2)
				if(c1 != [] and c2 != []):
					for x in range (0,len(c1)):
						Intersection_vector_1.append(c1[x])
						Intersection_vector_2.append(c2[x])	
		print 'Overlap Vector Length = ' + str(len(Intersection_vector_1))
		Corr_Data.append([-1 * i,dcor(Intersection_vector_1,Intersection_vector_2).find_correlation()])

		i = i + Shift_Step 
	
	#Used only for plotting correctly
	for s in LC_2:
		s.shift_segment(-1 * i)
	
	
	# Find Maximum correlation Correlation 
	print Corr_Data
	Shift = max(Corr_Data, key=operator.itemgetter(1))
	return Shift

'''
Returns a list of light curve segments fitting the points of 2 light curves
'''
def fit_curve(time, lc_a, lc_b, MAX_GAP):

# Perform segmentation
	T_Segments = []
	i = 0
	while (i < (len(time)-1)):
		temp1 = []
		temp2 = []
		tempt = []		
		while (i < (len(time)-1)) and (time[i+1]-time[i]<MAX_GAP):
			temp1.append(lc_a[i])
			temp2.append(lc_b[i])
			tempt.append(time[i])
			i = i + 1
		if(len(temp1)>5):	# Each segment must have at least 5 observational points	

		# s = smoothing parameter, if s = 0 curve passes through all the points.        //k = 3, s = 0.04
		# interpolate.UnivariateSpline(TimeScale, Intensity, k {Degree of Spline}, s {Smoothning Factor} )
			LC_1.append( Segment(interpolate.UnivariateSpline(tempt, temp1), temp1, tempt) ) 
			LC_2.append( Segment(interpolate.UnivariateSpline(tempt, temp2), temp2, tempt) ) 
		
			T_Segments.append(tempt)
		i = i + 1

	return LC_1, LC_2, T_Segments

if '__main__':
	
	# Load Data from file	
	time, lc_a, er_a, lc_b, er_b = load_data()

	# Introduce error by mutiplying given error with a standard normal distribution
	for i in range (0,len(lc_a)):
		lc_a[i] = lc_a[i] + random.gauss(0,1) * er_a[i]
		lc_b[i] = lc_b[i] + random.gauss(0,1) * er_b[i]
	
	# Compute mean value of both light curves and force it to zero	(not needed really)
	avg_a = sum(lc_a) / float(len(lc_a))
	avg_b = sum(lc_b) / float(len(lc_a))	

	for i in range (0,len(lc_a)):
		lc_a[i] = lc_a[i] - avg_a
		lc_b[i] = lc_b[i] - avg_b	


	# Perform Segmentation of Light Curve
	MAX_GAP = 60		# In Modified Heliocentric Julian Date	
	
	LC_1 = []
	LC_2 = []
	
	# Perform Segmentation and Spline fitting
	LC_1, LC_2, T_Segments = fit_curve(time, lc_a, lc_b, MAX_GAP)

	Max_Shift  = int((max(time)-min(time))/2)	
	Shift_Step = 2	

	# Find Distance Correlation 
	Shift = find_shift(LC_1, LC_2, Max_Shift, Shift_Step)

	print 'Maximum Observed Shift'
	print Shift

	print 'Shift : '+ str(Shift[0]) + ' mhjd days with correlation ' + str(Shift[1])
	

	#Plotting
	plot_data_points(time, lc_a, lc_b, avg_a, avg_b)
	for c1 in LC_1:
		plot_segment(c1, avg_a)
	for c2 in LC_2:
		plot_segment(c2, avg_b)
	show_plots()
