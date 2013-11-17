/* Exercise to convert a simple serial code for a brute force
   largest prime number search into OpenMP and MPI (32-bit, int version).

   All prime numbers can be expressed as 6*k-1 or 6*k+1, k being an
   integer. We provide the range of k to probe as macro parameters
   KMIN and KMAX (see below).

   Check the parallel code correctness - it should produce the same largest prime
   number as the serial version, for the same range KMIN...KMAX.

   Try to make the parallel code as efficient as possible.

   Your speedup should be close to the number of threads/ranks you are using.

OpenMP instructions:

   The code should print the number of threads used.

   Use "default(none)" in parallel region(s).

   Used OpenMP directives:
 - parallel
 - for schedule
 - critical
 - single

MPI instructions:

* Use simple static way of distributing the workload.
* All ranks (including rank 0) should participate in computations.
* The results should be correct for any number of ranks (so there should be no 
  assumption that the number of Monte Carlo steps is integer dividable by the
  number of ranks. As a result, some ranks might get a different amount of work.)
* Reduce the number of communications to a minimum.
* It should be the rank 0 responsibility to call gettimeofday() function, and print
  all the messages.
* Use the following functions:
 - MPI_Init
 - MPI_Finalize
 - MPI_Comm_rank
 - MPI_Comm_size
 - MPI_Reduce

Compiling instructions:

 - Serial code:
  icc -O2 primes.c -o primes

 - OpenMP code:
  icc -openmp -O2 primes_omp.c -o primes_omp


 - MPI code:
  mpicc -O2 primes_mpi.c -o primes_mpi
   To run:
  mpirun -np 24 ./primes_mpi

*/

#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// Range of k-numbers for primes search:
#define KMIN 100000000
// Should be smaller than 357,913,941 (because we are using signed int)
#define KMAX 101000000

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
  int devid, devcount, error, success;
  int xmax, ymax, x, y, k, j;

  gettimeofday (&tdr0, NULL);


  xmax = 0;
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
	      // Tpo be a success, the modulus should not be equal to zero:
	      success = x % y;
	      if (!success)
		break;
	    }

	  if (success && x > xmax)
	    {
	      xmax = x;
		  printf("xmax: %d\n", xmax);
	    }
	}
    }

  gettimeofday (&tdr1, NULL);
  tdr = tdr0;
  timeval_subtract (&restime, &tdr1, &tdr);
  printf ("maxval: %d\n", xmax);
  printf ("time: %e\n", restime);
  //--------------------------------------------------------------------------------



  return 0;

}
