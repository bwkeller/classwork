#!/usr/bin/python
#Solution for Part 2 of PHYS785 Problem Set 1
import numpy as np
import matplotlib.pyplot as plt

def lambdaH(temperature): #See Penston 1970 pp. 776
	return 1.634e-11*3.23e-10*np.sqrt(temperature)*(1+17500/temperature)*np.exp(-118000/temperature)
	return temperature
def lambdaCII(temperature, Cfrac=3.75e-4): #See Dalgarno & McCray 1972 Eq. 3.
	return Cfrac*7.9e-20/np.sqrt(temperature)*np.exp(-92/temperature)
def lambdaSII(temperature, Sfrac=1.4e-5): #See Penston 1970 pp. 780
	return Sfrac*8.44e-18/np.sqrt(temperature)*np.exp(-21400/temperature)
def lambdaTotal(temperature):
	return lambdaH(temperature)+lambdaCII(temperature)+lambdaSII(temperature)


if __name__ == "__main__":
	temperature = np.power(10, np.arange(1,6,0.01))
	plt.loglog(temperature, lambdaH(temperature), lw=2, label='Hydrogen')
	plt.loglog(temperature, lambdaCII(temperature), lw=2, label='Carbon II')
	plt.loglog(temperature, lambdaSII(temperature), lw=2, label='Sulfur II')
	plt.loglog(temperature, lambdaTotal(temperature), 'k--', lw=3, label='Total')
	plt.ylim(1e-28,1e-17)
	plt.legend(loc='best', title='Cooling Component')
	plt.xlabel("Temperature $(K)$")
	plt.ylabel(r"$\Lambda/n^2 (erg$ $cm^3 s^{-1})$")
	plt.show()
