
#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

// Range of k-numbers for primes search:
#define KMIN 100000000
// Should be smaller than 357,913,941 (because we are using signed int)
#define KMAX 101000000

// A quick macro for defining how much work should be done by each rank
#define RANGE ((KMAX-KMIN)/devcount)

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
  MPI_Init(&argc, &argv); //Initialize MPI
  struct timeval  tdr0, tdr1, tdr;
  double restime;
  int devid, devcount, error, success;
  int xmax, ymax, x, y, k, j, rootxmax;
  MPI_Comm_rank(MPI_COMM_WORLD, &devid); //Grab the current rank
  MPI_Comm_size(MPI_COMM_WORLD, &devcount); //Grab the total number of ranks
  int modulus = 0;

  //We only want rank 0 to do timings
  if (devid == 0)
  {
	  gettimeofday (&tdr0, NULL);
  }
  //Check to make sure that we aren't accidentally losing range 
  //because KMAX-KMIN isn't divisible by the number of ranks
  if (devid == devcount-1)
  {
	  modulus = (KMAX-KMIN) % devcount;
  }


  xmax = 0;
  //Test a range of numbers depending on your MPI rank
  for (k=KMIN+devid*RANGE; k<=KMIN+(devid+1)*RANGE+modulus; k++)
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
	    }
	}
    }

  //Use the builtin MPI_MAX to obtain the maximum of all the rank's local
  //xmax values and store the result in the root's variable rootxmax
  MPI_Reduce(&xmax, &rootxmax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  //Once again, leave timings and printings to the root rank.
  if(devid == 0)
  {
	  gettimeofday (&tdr1, NULL);
	  tdr = tdr0;
	  timeval_subtract (&restime, &tdr1, &tdr);
	  printf ("maxval: %d\n", rootxmax);
	  printf ("time: %e\n", restime);
  }
  //--------------------------------------------------------------------------------



  //Tear down MPI
  MPI_Finalize();
  return 0;

}
