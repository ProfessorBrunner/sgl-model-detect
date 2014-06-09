import numpy as np
import pylab as plt
import matplotlib
import random
import operator

#Probability of a seasonal effect being imprinted
Prob_Seasonal = 0.03

def C_factory(P, n=2, V_type="clamped"):
    """ Returns a b-spline curve C(t) configured with P, V_type and n.
    The knot vector will be created according to V_type
    
    Parameters
    ==========
    - P (list of D-tuples of reals) : List of de Boor points of dimension D.
    - n (int) : degree of the curve
    - V_type (str): name of the knit vector type to create.
    
    Returns
    =======
    A D-dimensionnal B-Spline Curve.
    """
    
    m = len(P)    # the number of points in P    
    D = len(P[0]) # the dimension of a point (2D, 3D)
        
    # Create the knot vector
    V = make_knot_vector(n, m, V_type)
    # TODO: check the validity of the input knot vector.
    # TODO: create an initial Vector Point.
    
    b_n = basis_factory(n)
    
    def S(t, d):
        """ The b-spline funtion, as defined in eq. 3. """
        out = 0.
        for i in range(m): #: Iterate over 0-indexed point indices
            out += P[i][d]*b_n(t, i, V)
        return out
    
    def C(t):
        """ The b-spline curve, as defined in eq. 4. """
        out = [0.]*D           #: For each t we return a list of D coordinates
        for d in range(D):     #: Iterate over 0-indexed dimension indices
            out[d] = S(t,d)
        return out
    
    C.P = P                   #: The control polygone
    C.V = V                   #: The knot vector used by the function
    C.spline = S              #: The spline function.
    C.basis = b_n             #: The highest degree basis function. Useful to do some plotting.
    C.min = V[0]              #: The domain of definition of the function, lower bound for t
    C.max = V[-1]             #: The domain of definition of the function, upper bound for t
    C.endpoint = C.max!=V[-1] #: Is the upper bound included in the domain.
    return C

def make_knot_vector(n, m, style="clamped"):
    """
    Create knot vectors for the requested vector type.
    
    Parameters
    ==========
    - n (int) : degree of the bspline curve that will use this knot vector
    - m (int) : number of vertices in the control polygone
    - style (str) : type of knot vector to output
    
    Returns
    =======
    - A knot vector (tuple)
    """
    if style != "clamped":
        raise NotImplementedError
        
    total_knots = m+n+2                           
    outer_knots = n+1                             
    inner_knots = total_knots - 2*(outer_knots)   
    # Now we translate eq. 5:
    knots  = [0]*(outer_knots)
    knots += [i for i in range(1, inner_knots)]
    knots += [inner_knots]*(outer_knots)
    
    return tuple(knots) 

def basis_factory(degree):
    """ Returns a basis_function for the given degree """
    if degree == 0:
        def basis_function(t, i, knots):
            """The basis function for degree = 0 as per eq. 7"""
            t_this = knots[i]
            t_next = knots[i+1]
            out = 1. if (t>=t_this and t<t_next) else 0.         
            return out
        
    else:
        def basis_function(t, i, knots):
            """The basis function for degree > 0 as per eq. 8"""
            out = 0.
            t_this = knots[i]
            t_next = knots[i+1]
            t_precog  = knots[i+degree]
            t_horizon = knots[i+degree+1]            

            top = (t-t_this)
            bottom = (t_precog-t_this)
     
            if bottom != 0:
                out  = top/bottom * basis_factory(degree-1)(t, i, knots)
                
            top = (t_horizon-t)
            bottom = (t_horizon-t_next)
            if bottom != 0:
                out += top/bottom * basis_factory(degree-1)(t, i+1, knots)
         
            return out
        
    # Add information #
    basis_function.lower = None if degree==0 else basis_factory(degree-1)
    basis_function.degree = degree
    return basis_function


def generateFile(filename, D ):

	#Write to file
	output = open(filename, 'w')
	output.write('mhjd\tmag_A\tmagerr_A\tmag_B\tmagerr_B\ttelescope\n')
	output.write('====\t=====\t========\t=====\t========\t=========\n')
	
	for i in range (0,len(D)):
		output.write(str(D[i][0])+'\t')
		output.write(str(D[i][1])+'\t')
		output.write(str(D[i][2])+'\t')
		output.write(str(D[i][3])+'\t')
		output.write(str(D[i][4])+'\t')
		output.write('None'+'\n')

	output.close()

def trimCurve(Flux1, error1, time1, Flux2, error2, time2):
	X=[]
	while (time1[0] < time2[0]):		
		del Flux1[0]
		del time1[0]
		del error1[0]
	while (time2[0] < time1[0]):		
		del Flux2[0]
		del time2[0]
		del error2[0]

	for i in range (0,min(len(Flux1),len(Flux2))):
		X.append([time2[i],Flux1[i],error1[i],Flux2[i], error2[i]])
	return X
		

def getError(Mean,Variance,tsteps): 
	error = []		
	for i in range(0, tsteps):
		error.append(random.gauss(Mean,Variance))
	return error	


'''
Plots the bezier curve
'''
def plot_data_points(time, lc, C):
	plt.plot(time, lc, marker="o")
	plt.plot(*zip(*C.P), alpha=0.3)

def show_plots():
	plt.legend(loc=4)
	plt.xlabel('Modified Heliocentric Julian Date (MHJD)')
	plt.ylabel('Magnitude (relative)')
	plt.show()



if '__main__':

	'''        
	t1 = input('Enter curve start time: ')
	t2 = input('Enter curve end time: ')

	i1 = input('Enter curve Max Intensity: ')
	i2 = input('Enter curve Min Intensity: ')
	'''

	plt.axis([1000, 1200, 8, 10])
	#num = input('Please enter number of points: ')
	pts = plt.ginput(4) # it will wait for three clicks


	LC1 = C_factory(P=pts, n=2, V_type="clamped")
	lenPts1 = int(max(pts, key=operator.itemgetter(1))[0])- int(min(pts, key=operator.itemgetter(1))[0])
	time = [t for t in np.linspace(int(LC1.min), int(LC1.max), lenPts1, endpoint=int(LC1.endpoint))]
	curvepts = [ LC1(s) for s in time ]
	Flux1 = []
	time1 = []
	for i in range(0,len(curvepts)):
		Flux1.append(curvepts[i][1])
		time1.append(curvepts[i][0])	       
	error1 = getError(1, 1, lenPts1)
			   
	plot_data_points(time1,Flux1,LC1)


	#num = input('Please enter number of points: ')
	pts = plt.ginput(4) # it will wait for three clicks

	LC2 = C_factory(P=pts, n=2, V_type="clamped")
	lenPts2 = int(max(pts, key=operator.itemgetter(1))[0])- int(min(pts, key=operator.itemgetter(1))[0])
	time = [t for t in np.linspace(int(LC2.min), int(LC2.max), lenPts2, endpoint=int(LC2.endpoint))]

	curvepts = [ LC2(s) for s in time ]
	Flux2 = []
	time2 = []
	for i in range(0,len(curvepts)):
		Flux2.append(curvepts[i][1])
		time2.append(curvepts[i][0])	       

	error2 = getError(2, 1, lenPts2)

	plot_data_points(time2,Flux2,LC2)
	
	print 'Time1'
	print time1
	print 'Time2'
	print time2
	
	X = trimCurve(Flux1, error1, time1, Flux2, error2, time2 )
	generateFile('Bezier_Curve.rdb',X)
	show_plots()
