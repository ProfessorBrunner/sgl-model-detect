import numpy as np
from dcor import * 

import pylab


def stochasticIntegral(tau):
    steps = 5 # Default number of sampling points for integral evaluation
    
    # Normal parameters
    mean = 0.0
    sigma = 1.0
    
    # Our integral sample points
    sample = np.arange(steps)
    
    # Our integration is simply summing up the random dt and the function value

    values = np.exp(-sample/tau) * np.random.normal(mean, sigma, steps)

    return (np.sum(values))

def initTimescales(mass, radius, alpha):
    # Renormalize to base units
    mass /= 1.0E8
    radius /= 100.0

    # The B. Kelly model has three timescales
    tau = np.zeros(3)
    
    # First we compute the light-crossing timescale
    tau[0] = 1.1 * mass * radius
    
    # Second, we compute the orbital timescale
    tau[1] = 104.0 * mass * (radius**1.5)
    
    # Third we calculate the thermal timescale
    diny = 365 # days in a year.
    tau[2] = 4.6 * (0.01/alpha) * mass * (radius**1.5) * diny

    return tau

def generateLC(mass, radius, alpha, b, sigma, tsteps):

    # Define parameters:
    ntimes = 3 # Number of timescales
    
    # First get timescales
    tau = initTimescales(mass, radius, alpha)
    
    # Initialize arrays
    fluxes  = np.zeros((ntimes,tsteps))
    ssi = np.zeros(ntimes)
    
    # Initial random flux values
    fluxes[:,0] = np.random.normal(b * tau, sigma * (tau/2.0)**0.5)

    # Common calculations    
    bt = b * tau
    exp_tau = np.exp(-5.0/tau)

    for i in range(1, tsteps):
        # Do our three stochastic integrals
        for j in range(ntimes):
            ssi[j] = sigma * stochasticIntegral(tau[j])
            
        fluxes[:,i] = exp_tau * fluxes[:, i - 1] + bt * (1.0 - exp_tau) + ssi
 
    return fluxes


tsteps = 512
mass = 1.0E8
radius = 100
sigma = 1.0
b = 0.
alpha = 0.01

flux = generateLC(mass, radius, alpha, b, sigma, tsteps)
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

L1 = []						# Data points for Light Curve
L2 = []						# Data points for Shifted Curve
Itr= [] 
DCor = []
L1_tstart = 200
Curve_Shift = 500

L1_tend   = tsteps + L1_tstart -1
L2_tstart = L1_tstart +Curve_Shift 
L2_tend   = tsteps + L2_tstart -1

day = 0
for i in range(0,2000):
	if(day<L1_tstart):	
		L1.append(0)
	elif (day>L1_tend):
		L1.append(0)
	else:
		L1.append(flux[0,:][day-L1_tstart])
	
	if(day<L2_tstart):
		L2.append(0)
	elif (day>L2_tend):
		L2.append(0)
	else:
		L2.append(flux[0,:][day-L2_tstart])
	Itr.append(day)
	day+=1


f = open('Output.txt', 'w')
Tau = []
for i in range (L1_tstart, L1_tstart + 1000):
	Sample = []
	for x in range(i, i + tsteps):
		Sample.append(L2[x])	
	X = dcor(flux[0,:],Sample)
	DCor.append(X.find_correlation())
	Tau.append(i - L1_tstart)
	f.write(str(X.find_correlation())+ ' '+ str(i-L1_tstart) +'\n')
	print i-L1_tstart
f.close()


pylab.plot(Tau, DCor, label='Distance Correlation')
pylab.legend()
pylab.xlabel('Shift (Tau)')
pylab.ylabel('Distance Correlation')
pylab.title('Light Curve Shifting')

#pylab.plot(Itr, L1, label='Light Curve 1')
#pylab.plot(Itr, L2, label='Light Curve 2')
#pylab.legend()
#pylab.xlabel('Day')
#pylab.ylabel('Y = f(x)')
#pylab.title('Light Curve Shifting')

pylab.show()
