#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

def planckfunction(wl, temp):
	return 1.19e-16/np.power(wl, 5.0)*(1.0/(np.exp((0.0144)/(wl*temp))-1.0))

if __name__ == "__main__":
	Irange = np.concatenate((np.arange(1e-9, 1e-6, 1e-9), np.arange(1e-6, 1e-3, 1e-6), np.arange(1e-3, 1, 1e-3), np.arange(1, 1000, 1)))
	I30000 = np.vectorize(lambda x: planckfunction(x, 30000))
	I6000 = np.vectorize(lambda x: planckfunction(x, 6000))
	I2500 = np.vectorize(lambda x: planckfunction(x, 2500))
	plt.loglog(Irange, I2500(Irange))
	plt.loglog(Irange, I6000(Irange))
	plt.loglog(Irange, I30000(Irange))
	plt.ylim(1,1e14)
	plt.axvline(3.5e-7)
	plt.axvline(8.5e-7)
	plt.show()
