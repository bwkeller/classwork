#!/usr/bin/python
# Ben Keller
# Physics 746
# Assignment 2, Exercise 1 Solution
import numpy as np
import matplotlib.pyplot as plt

#Physical Constants
h = 6.63e-34
c = 3e8
kB = 1.38e-23

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

#Simple Black Body Spectrum
def planck(T, nu):
	return 2*h*np.power(nu, 3)/(c*c*(np.exp(h*nu/(kB*T)) - 1 ))

#The temperature function for a grey, plane-parallel atmosphere
def temperature(T_eff, tau):
	return np.power((0.75*np.power(T_eff, 4)*(tau+2./3)), 0.25)

#The first integral (of dx)
def integrand_x(tau, x):
	return np.power(x, -2)*np.exp(-1*x*tau)

#The second integral (of dtau)
def integrand_tau(T_eff, tau, nu, end):
	return planck(temperature(T_eff, tau), nu)*simpson_integrate(lambda x: integrand_x(tau, x), 1, end, 1)

#Full Spectrum for our grey, plane-parallel atmosphere
def grey_flux(T_eff, nu, end):
	return 2*np.pi*simpson_integrate(lambda x: 
	integrand_tau(T_eff, x, nu, end), 0, end, 1)

if __name__ == "__main__":
	freq = np.power(10, np.arange(10, 16, 0.1))
	plt.loglog(freq, np.pi*planck(1000, freq), 'r-', label=r'$T_{eff}$ = 1000K')
	plt.loglog(freq, grey_flux(1000, freq, 1e3), 'r--')
	plt.loglog(freq, np.pi*planck(5000, freq), 'g-', label=r'$T_{eff}$ = 5000K')
	plt.loglog(freq, grey_flux(5000, freq, 1e3), 'g--')
	plt.loglog(freq, np.pi*planck(10000, freq), 'b-', 
	label=r'$T_{eff}$ = 10000K')
	plt.loglog(freq, grey_flux(10000, freq, 1e3), 'b--')
	plt.legend(loc='best')
	plt.xlabel("Frequency, log10 scaling")
	plt.ylabel("Flux")
	plt.show()
