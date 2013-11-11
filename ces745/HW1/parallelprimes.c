// CES 745 Homework #1 Ben Keller 1178881
#include <stdio.h>
#include <omp.h>
#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// Range of k-numbers for primes search:
#ifndef KMIN
#define KMIN 100000000
#endif
// Should be smaller than 357,913,941 (because we are using signed int)
#ifndef KMAX
#define KMAX 101000000
#endif

/* Subtract the `struct timeval' values X and Y,
   storing the result in RESULT.
   Return 1 if the difference is negative, otherwise 0.  */

// It messes up with y!

	int
timeval_subtract (double *result, struct timeval *x, struct timeval *y)
{
	struct timeval result0;

	/* Perform the carry for the later subtraction by updating y. */
	if (x->tv_usec < y->tv_usec) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_usec - y->tv_usec > 1000000) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}

	/* Compute the time remaining to wait.
	   tv_usec is certainly positive. */
	result0.tv_sec = x->tv_sec - y->tv_sec;
	result0.tv_usec = x->tv_usec - y->tv_usec;
	*result = ((double)result0.tv_usec)/1e6 + (double)result0.tv_sec;

	/* Return 1 if result is negative. */
	return x->tv_sec < y->tv_sec;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


int main (int argc,char **argv)
{
	struct timeval  tdr0, tdr1, tdr;
	double restime;
	int devcount, error, success;
	int xmax, ymax, x, y, k, j;
	int devxmax; //This is a thread-private version of xmax (which is the global max)

	gettimeofday (&tdr0, NULL);


	xmax = 0;
#pragma omp parallel default(none) shared(xmax, devcount) private(k, j, x, ymax, y, success, devxmax)
	{
		devcount = omp_get_num_threads(); //Current thread ID stored in devid
		devxmax = 0;//Zero out the thread-specific maximum
#pragma omp for schedule(dynamic)
		for (k=KMIN; k<=KMAX; k++)
		{
			// testing "-1" and "+1" cases:
			for (j=-1; j<2; j=j+2)
			{
				// Prime candidate:
				x = 6*k + j;
				// We should be dividing by numbers up to sqrt(x):
				ymax = (int)ceil(sqrt((double)x));

				// Primality test:
				for (y=3; y<=ymax; y=y+2)
				{
					// To be a success, the modulo should not be equal to zero:
					success = x % y;
					if (!success)
						break;
				}

				if (success && x > devxmax)
				{
					devxmax = x;
				}
			}
		}
		/* The comparison to check for the largest prime mustn't be done by 
		 * multiple threads simultaneously, or a race condition may occur */
#pragma omp critical
		{
			if (devxmax > xmax)
			{
				xmax = devxmax;
			}
		}
	}
/*This must only be run by a single thread, since we don't want a pile of 
 * printfs each with slightly different times */
#pragma omp single 
	{
		gettimeofday (&tdr1, NULL);
		tdr = tdr0;
		timeval_subtract (&restime, &tdr1, &tdr);
		printf ("threads: %d\n", devcount);
		printf ("maxval: %d\n", xmax);
		printf ("time: %e\n", restime);
	}
	//--------------------------------------------------------------------------------



	return 0;

}
