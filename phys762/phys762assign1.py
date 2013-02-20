#!/usr/bin/python
# Phsyics 762 Assignment 1 Question 3
# Ben Keller
# 1178881
# November 2 2011
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
	SmNd = [0.1847, 0.1963, 0.1980, 0.2061, 0.2715, 0.2879]
	NdNd = [0.511721, 0.511998, 0.512035, 0.512238, 0.513788, 0.514154]
	NdNdError = [0.000018, 0.000016, 0.000021, 0.000017, 0.000015, 0.000017]
	plt.errorbar(SmNd, NdNd, yerr=NdNdError, color='k', linestyle='--')
	plt.xlabel("Sm-147/Nd-144 Ratio")
	plt.ylabel("Nd-143/Nd-144 Ratio")
	plt.title("Samarium-Neodymium Isotope Dating of Lunar Basalt")
	print np.polyfit(SmNd, NdNd, 1)
	plt.savefig("phys762assign1.pdf")
	plt.show()
