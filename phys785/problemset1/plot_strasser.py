#!/usr/bin/python
#Generate the Strasser & Taylor 2004 style plots for part 1 of PHYS785 problem set 1.
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

if __name__ == "__main__":
	datfile = np.genfromtxt(argv[1], skiprows=1)
	LSR = datfile[:,0]
	Toff = datfile[:,1]
	Ton = datfile[:,2]
	exptau = datfile[:,3]
	dToff = datfile[:,4]
	dTon = datfile[:,5]
	dexptau = datfile[:,6]
	Ts = Toff/(1-exptau)
	dTs = dToff/(1-exptau)+dexptau/np.power(exptau-1, 2)
	SN = np.abs(Toff-Ton)/np.sqrt(dToff**2+dTon**2)
	highSN = np.where(SN > 3)[0]
	plt.subplot(311)
	plt.plot(LSR, Toff, 'k-')
	plt.title("$l = %s \degree, b = %s \degree$" % (argv[1].split('_')[0], argv[1].split('_')[1][:-4]))
	plt.ylabel(r'$T_{boff} (K)$')
	plt.xlim(-160,50)
	plt.subplot(312)
	plt.errorbar(LSR, exptau, yerr=dexptau, color='k', marker='.')
	plt.ylim(0,1.3)
	plt.xlim(-160,50)
	plt.ylabel(r'$<e^\tau>$')
	plt.subplot(313)
	plt.errorbar(LSR[highSN], Ts[highSN], color='k', linestyle='none', marker='.', yerr=dTs[highSN])
	plt.xlim(-160,50)
	plt.ylabel(r'$<T_S> (K)$')
	plt.xlabel(r'Velocity $(km/s)$')
	plt.tight_layout()
	plt.show()
