#!/usr/bin/python
# Physics 762 Assignment 2 Question 1
# Ben Keller
# 1178881
# November 8 2011
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
	V = np.arange(4,20)
	Am = [-2.31, -1.71, -1.11, -0.51, 0.09, 0.69, 1.29, 1.89, 2.24, 2.59, 2.94,
	3.29, 3.89, 4.49, 5.09, 5.69]
	plt.plot(V, Am, 'ko')
	(slope, int1) = np.polyfit(V[:8], Am[:8], 1)
	(slope, int2) = np.polyfit(V[-5:], Am[-5:], 1)
	plt.plot(V, V*slope+int1, 'r-', label='Fitted Slope: %1.1f Fitted Intercept: %1.2f' % (slope, int1))
	plt.plot(V, V*slope+int2, 'b-', label='Fitted Slope: %1.1f Fitted Intercept: %1.2f' % (slope, int2))
	plt.legend(loc=0)
	plt.ylabel(r'$\log_{10}A_M$')
	plt.xlabel(r'$V$')
	plt.title('Gas Cloud Extinction')
	plt.savefig('phys762assign2q1.pdf')
	plt.show()
