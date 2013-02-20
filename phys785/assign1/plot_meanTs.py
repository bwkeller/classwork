#!/usr/bin/python
#Generate the mean spin temperatures for problem 1 of problem set 1 for PHYS785
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

if __name__ == "__main__":
	datfile = np.genfromtxt('110.600_4.6946.dat', skiprows=1)
	LSR = datfile[:,0]
	Toff = datfile[:,1]
	Ton = datfile[:,2]
	exptau = datfile[:,3]
	dToff = datfile[:,4]
	dTon = datfile[:,5]
	dexptau = datfile[:,6]
	Ts1 = Toff/(1-exptau)
	dTs1 = dToff/(1-exptau)+dexptau/np.power(exptau-1, 2)
	SN1 = np.abs(Toff-Ton)/np.sqrt(dToff**2+dTon**2)
	highSN1 = np.where(SN1 > 3)[0]
	datfile = np.genfromtxt('111.612_0.3739.dat', skiprows=1)
	LSR = datfile[:,0]
	Toff = datfile[:,1]
	Ton = datfile[:,2]
	exptau = datfile[:,3]
	dToff = datfile[:,4]
	dTon = datfile[:,5]
	dexptau = datfile[:,6]
	Ts2 = Toff/(1-exptau)
	dTs2 = dToff/(1-exptau)+dexptau/np.power(exptau-1, 2)
	SN2 = np.abs(Toff-Ton)/np.sqrt(dToff**2+dTon**2)
	highSN2 = np.where(SN2 > 3)[0]
	datfile = np.genfromtxt('111.649_2.4038.dat', skiprows=1)
	LSR = datfile[:,0]
	Toff = datfile[:,1]
	Ton = datfile[:,2]
	exptau = datfile[:,3]
	dToff = datfile[:,4]
	dTon = datfile[:,5]
	dexptau = datfile[:,6]
	Ts3 = Toff/(1-exptau)
	dTs3 = dToff/(1-exptau)+dexptau/np.power(exptau-1, 2)
	SN3 = np.abs(Toff-Ton)/np.sqrt(dToff**2+dTon**2)
	highSN3 = np.where(SN3 > 3)[0]
	highSN = np.intersect1d(highSN1, np.intersect1d(highSN2, highSN3))
	Tsbar = Ts1 + Ts2 + Ts3
	dTsbar = np.sqrt((dTs1**2+dTs2**2+dTs3**2)/3)
	plt.errorbar(LSR[highSN], Tsbar[highSN], color='k', linestyle='none', marker='.', yerr=dTsbar[highSN])
	plt.ylabel(r'$<\bar T_S> (K)$')
	plt.xlabel(r'Velocity $(km/s)$')
	plt.tight_layout()
	plt.show()
