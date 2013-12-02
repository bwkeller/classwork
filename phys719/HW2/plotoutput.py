import matplotlib.pyplot as plt
import numpy as np

ncells=400
xmax=2
dx=2.0*xmax/ncells

if __name__ == "__main__":
	infile = np.genfromtxt("output.dat")
	plt.subplot(311)
	plt.plot(np.arange(-xmax,xmax,dx), infile[:,0])
	plt.xlabel('x')
	plt.ylabel('v')
	plt.subplot(312)
	plt.plot(np.arange(-xmax,xmax,dx), infile[:,1])
	plt.xlabel('x')
	plt.ylabel('p')
	plt.subplot(313)
	plt.plot(np.arange(-xmax,xmax,dx), infile[:,2])
	plt.xlabel('x')
	plt.ylabel(r'$\rho$')
	fig = plt.gcf()
	fig.set_size_inches(8.5,11)
	plt.savefig("output.pdf")