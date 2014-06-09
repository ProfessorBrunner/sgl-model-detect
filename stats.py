import sys
import random

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


if '__main__':
	'''	
	time, lc_a, er_a, lc_b, er_b = load_data()
	
	avg_a = sum(er_a)/len(er_a)
	avg_b = sum(er_b)/len(er_b)


	sum_a, sum_b = 0 , 0
	for i in range (0,len(er_a)):
		sum_a = sum_a + (avg_a - er_a[i])**2
		sum_b = sum_b + (avg_b - er_b[i])**2
	
	sum_a = sum_a / len(er_a)**2
	sum_b = sum_b / len(er_b)**2

	sum_a = sum_a**0.5
	sum_b = sum_b**0.5
	'''
	desired_mean = 2.
	print random.expovariate(1/desired_mean)


