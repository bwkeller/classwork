#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

#The integral differential element of the flux off the slab
def dslabflux(z, theta, tau0):
	return np.exp(-1.0*z*tau0/np.cos(theta))*np.cos(theta)*np.sin(theta)

#This function calculates the integral of a function evaluated from minimum to
#maximum using Simpson's rule, with a grid size of step
def simpson_integrate(function, minimum, maximum, step):
	integral = 0.0
	n = 1
	while maximum > minimum+n*step:
		if n % 2 == 0:
			integral += 2.0*function(minimum+n*step)
		else:
			integral += 4.0*function(minimum+n*step)
		n += 1
	return step/3.0*(function(minimum) + integral)

#The total integrated flux coming off the slab
def slabflux(z, tau0):
	return 2*np.pi*(0.5-simpson_integrate(lambda x: dslabflux(z, x, tau0), 0, 
	np.pi*0.5, 1e-4))-2*np.pi*(0.5-simpson_integrate(lambda x: dslabflux(1-z, 
	x, tau0), 0, np.pi*0.5, 1e-4))

if __name__ == "__main__":
	zrange = np.power(10, (np.arange(-3, 0, 0.001))) #z range [0.001, 1]
	piline = np.pi*np.ones(zrange.size)
	tau = (0.01, 0.1, 1, 10, 100) #tau0 range
	plt.plot(zrange , slabflux(zrange, tau[0]), 
	label=r"$\tau_0=%3.2f$" % (tau[0]))
	plt.plot(zrange , slabflux(zrange, tau[1]),
	label=r"$\tau_0=%3.2f$" % (tau[1]))
	plt.plot(zrange , slabflux(zrange, tau[2]),
	label=r"$\tau_0=%d$" % (tau[2]))
	plt.plot(zrange , slabflux(zrange, tau[3]),
	label=r"$\tau_0=%d$" % (tau[3]))
	plt.plot(zrange , slabflux(zrange, tau[4]),
	label=r"$\tau_0=%d$" % (tau[4]))
	plt.xlabel(r"Slab Depth")
	plt.ylabel(r"Flux$/B_\nu$")
	plt.title("Exercise 2.4")
	plt.legend(loc=2)
	plt.show()
