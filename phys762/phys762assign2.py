#!/usr/bin/python
# Physics 762 Assignment 2 Question 5
# Ben Keller
# 1178881
# November 8 2011

import numpy as np
import matplotlib.pyplot as plt

M0 = 7.36e36
k = 6.2e19
G = 6.67e-11
Msun = 1.989e30
parsec = 3.086e16

def mass(radius):
	return k*radius+M0

def velocity(radius):
	return np.sqrt(G*(k+M0/radius))

if __name__ == "__main__":
	radialrange = parsec*np.arange(0.01, 1000, 0.01)
	massrange = mass(radialrange)
	velocityrange = velocity(radialrange)
	plt.loglog(radialrange/parsec, massrange/Msun, 'k-')
	plt.xlabel(r'$\log_{10}r (pc)$')
	plt.ylabel(r'$\log_{10}M_r (M\odot)$')
	plt.title('Mass Distribution for Inner Galaxy')
	plt.savefig('phys762assign2d.pdf')
	plt.show()
	plt.cla()
	plt.loglog(radialrange/parsec, velocityrange/1000.0, 'k-')
	plt.xlabel(r'$\log_{10}r (pc)$')
	plt.ylabel(r'$\log_{10}v (kms^{-1})$')
	plt.title('Velocity Distribution for Inner Galaxy')
	plt.savefig('phys762assign2e.pdf')
	plt.show()
