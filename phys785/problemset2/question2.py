#!/usr/bin/python
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

if __name__ == "__main__":
	masses = np.genfromtxt("20Kmasses.dat")
	#Apply the low-mass threshold (0.4 Msol is when incompleteness becomes
	#significant, but Johnstone et al. 2000 found a change in slope at 0.6 Msol
	masses = masses[np.where(masses > 0.6)]
	#Take a mass histogram to find dN/dM.  6 bins is the maximum number that has no
	#empty bins
	(count, bins) = np.histogram(np.log10(masses), 6)
	binerr = (bins[1:]-bins[:-1])/2
	bincenters = bins[:-1]+binerr
	#Relative Uncertainty in N is sqrt(N)/N
	counterr = 1/np.sqrt(count)
	count = count.astype(np.float64)/np.power(10, bincenters)
	#Generate a power-law fit using a uncertainty-weighted least-squares fit of
	#log(N) to log(M)
	fitfunc = lambda p, x: p[0]+p[1]*x
	errfunc = lambda p, x,y, yerr,xerr: (y-fitfunc(p,x))/np.sqrt(xerr**2+yerr**2)
	fit = leastsq(errfunc, [-0.2,6], args=(bins[:-1], np.log10(count), counterr, binerr),
			full_output=1)
	(intercept, slope) = fit[0]
	#Print the slope and the covariance matrix's inner product with the slope (uncertainty)
	print "Slope of Power Law Fit %3.2f +/- %3.2f" % (slope, np.sqrt(fit[1][0][0]))
	plt.plot(bincenters, bincenters*slope+intercept, 'k--', label='Power-Law Fit')
	plt.errorbar(bincenters, np.log10(count), linestyle='none', color='red', marker='o', 
			yerr=counterr, xerr=binerr, label='Clump Observation')
	plt.xlabel('$log_{10}$ Mass $(M_\odot)$')
	plt.ylabel(r'$log_{10}$ $dN/dM$ $(M^{-1}_\odot)$')
	plt.legend()
	plt.show()
