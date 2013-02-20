#!/usr/bin/python
# Ben Keller
# October 2011
# Solution for Question 2 in Physics 761 Assignment 2
import numpy
import matplotlib.pyplot as plt

#The Lane Emden Equation, reformulated using the subsitution z=Dn' to be
#z' = -Dn^n - 2z/x
def LaneEmden(xi, Dn, z, n):
	return -1.0*numpy.power(Dn, n)-(2.0/xi)*z


#This function is simply the Runge-Kutta 4th order method for solving ODEs.
#It solves for y' = func(x, y) at the point x=xf, and returns all the 
#intermediate values
def rk4(func, h, x0, y0, z0, xf):
	yn = y0
	xn = x0
	zn = z0
	yvals = []
	xvals = []
	while xn < xf:
		k1 = h*func(xn, yn, zn)
		k2 = h*func(xn+0.5*h, yn+0.5*k1, zn)
		k3 = h*func(xn+0.5*h, yn+0.5*k2, zn)
		k4 = h*func(xn+h, yn+k3, zn)
		zn += (k1 + k2 + k3 + k4)/6.0
		yn += h*zn
		xn += h
		yvals.append(yn)
		xvals.append(xn)
	return (xvals, yvals)

if __name__ == "__main__":
	#Define the n=1.5 and n=3.0 Lane Emden equations
	lane1 = lambda x, y, z: LaneEmden(x, y, z, 1.5)
	lane2 = lambda x, y, z: LaneEmden(x, y, z, 3.0)
	#Calculate the points for the two solutions
	(x1, y1) = rk4(lane1, 1e-5, 1-((1e-10)/6), 1, 0, 5)
	(x2, y2) = rk4(lane2, 1e-5, 1-((1e-10)/6), 1, 0, 6)
	#Plot the results
	plt.plot(x1, y1, 'k-', label="n=1.5")
	plt.plot(x2, y2, 'k--', label="n=3")
	plt.legend()
	plt.ylabel("$D_n$")
	plt.xlabel("$\\xi$")
	plt.title("Solutions to the Lane-Emden Equation")
	plt.ylim((0,1))
	plt.savefig("phys761assign2.pdf")
	plt.show()
