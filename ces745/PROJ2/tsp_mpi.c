/* Exercise to solve the Traveling Salesman Problem (TSP) using
 a brute force Monte Carlo approach. Your task is to convert the serial
 computation to a parallel computation. The code has to print the
 current shortest total traveling distance found fairly frequently, so 
 you can monitor its progress.

 We use random number generator rand_r() specifically because it is thread-safe
 (it doesn't keep anything in a global state; everything is defined by the 
 variable "seed", which can be made private to the thread). 

 MPI: let rank 0 compute the initial distance matrix, and send it to other ranks.
 The timing (both serial and MPI code) should be done from right before sending 
 the distance matrix to other ranks, until right before printing the result by
 rank 0 (so it should include all MPI operations on rank 0, except for MPI_Init and 
 MPI_Finalize). The rank 0 has to print both the shortest distance and the
 corresponding itinerary.

 MPI: You should get the correct result for any N_MC (that is, N_MC doesn't have
 to be integer dividable by the number of ranks.

 Make sure the parallel code produces correct results, both in terms of
 the smallest distance found and the corresponding itinerary. To do
 this, run both the serial and parallel versions of the code for a small
 number of cities, N_CITIES - say 10 or 11, and identical total number
 of Monte Carlo steps, N_MC. Make sure N_MC is a few times larger than
 the factorial of (N_CITIES-1), so despite the random nature of the
 search you are almost guaranteed to find the exact solution.

 Your speedup should be close to the number of threads/ranks you are using.
 Say, on orca development nodes (24 threads) speedup should be >~23.

 For convenience, you can develop the code with smaller numbers (N_CITIES=10,
 N_MC=1e7). For final timing tests use N_CITIES=11 and N_MC=1e8.


 - Don't change the way the cities coordinates are initially randomly
 generated on the host (this will make it difficult to judge
 the code performance).

 - Don't forget that int can only handle integers up to ~4 billion. 

 - As this is Monte Carlo simulation, don't forget that in a parallel
 version each thread/rank has to start with a unique seed number.

 - MPI:
mpicc -O2 tsp_mpi.c -o tsp_mpi


To run (on orca devel nodes):
 - MPI:
mpirun -np 24 ./tsp_mpi


*/

#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

// Number of cities:
#define N_CITIES 11

// Maximum number of Monte Carlo steps:
#define N_MC (long long int) 1e8


// Size of the square region (km) where the cities are randomly scattered:
#define SIZE 1000

// Cities coordinates:
float x[N_CITIES], y[N_CITIES];
// Distance array:
float dist[N_CITIES][N_CITIES];


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
  int devid, devcount, error, i, j, k, l, l1, ltemp;
  int perm[N_CITIES], perm_min[N_CITIES];
  float d;
  MPI_Comm_rank(MPI_COMM_WORLD, &devid); //Grab the current rank
  MPI_Comm_size(MPI_COMM_WORLD, &devcount); //Grab the total number of ranks
  int modulus = 0;

  //We only want rank 0 to generate the distance matrix & timing
  if (devid == 0)
  {
	  // Generating randomly scattered cities
	  for (i=0; i<N_CITIES; i++)
		{
		  x[i] = (float)SIZE * (float)rand()/((float)RAND_MAX+1.0);
		  y[i] = (float)SIZE * (float)rand()/((float)RAND_MAX+1.0);
		}

	  // Computing the distance matrix:
	  for (i=N_CITIES-1; i>=0; i--)
		{
		  for (j=0; j<N_CITIES; j++)
		{
		  if (j < i)
			dist[i][j] = sqrt(pow(x[j]-x[i],2) + pow(y[j]-y[i],2));
		  else if (j == i)
			dist[i][j] = 0.0;
		  else
			dist[i][j] = dist[j][i];
		}
		}
	  
	  gettimeofday (&tdr0, NULL);
  }
  MPI_Bcast(x, N_CITIES, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(y, N_CITIES, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dist, N_CITIES*N_CITIES, MPI_FLOAT, 0, MPI_COMM_WORLD);


  //Check to make sure that we aren't accidentally losing range 
  //because KMAX-KMIN isn't divisible by the number of ranks
  if (devid == devcount-1)
  {
	  modulus = N_MC % devcount;
  }

  float d_min = 1e30;
  float rootd_min; 
  unsigned int seed = 111+devid;

  // We always start from the same (0-th) city:
  perm[0] = 0;

  // Cycle for Monte Carlo steps:
  for (k=0; k<N_MC/devcount + modulus; k++)
    {
      // Generating a random permutation:

      // Initially we have an ordered list:
      for (l=1; l<N_CITIES; l++)
	perm[l] = l;

      // Then we reshuffle it randomly, starting with l=1:
      d = 0.0;
      for (l=1; l<N_CITIES-1; l++)
	{
	  // This generates a random integer in the range [l ... N_CITIES-1]:
	  l1 = l + rand_r(&seed) % (N_CITIES-l);

	  // Swapping the l and l1 cities:
	  ltemp = perm[l];
	  perm[l] = perm[l1];
	  perm[l1] = ltemp;

	  // At this point, cities in perm[l-1] and perm[l] have already been reshuffled, so we 
	  // can compute their contribution to the total distance:
	  d = d + dist[perm[l-1]][perm[l]];
	}

      // At the final leg we are coming back to the original (0-th) city:
      d = d + dist[perm[N_CITIES-1]][0];
      //!!! It was a bug:
      d = d + dist[perm[N_CITIES-1]][perm[N_CITIES-2]];

      // Finding globally shortest TSP distance:
      if (d < d_min)
	{
	  d_min = d;
	  printf ("%d %f ", k, d_min);
	  // Printing the itinerary corresponding to the current smallest distance:
	  for (l=0; l<N_CITIES; l++)
	    {
	      printf ("%d ", perm[l]);
	      // Memorizing the shortest itinernary:
	      perm_min[l] = perm[l];
	    }
	  printf ("\n");
	}

    }  
  
  //Use the builtin MPI_MAX to obtain the maximum of all the rank's local
  //xmax values and store the result in the root's variable rootxmax
  MPI_Allreduce(&d_min , &rootd_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  //Once again, leave timings and printings to the root rank.
  if(devid == 0)
  {
	  gettimeofday (&tdr1, NULL);
	  tdr = tdr0;
	  timeval_subtract (&restime, &tdr1, &tdr);
	  printf ("Shortest total distance: %f\n", d_min);

	  // Timing for the region to be parallelized:
	  printf ("%e\n", restime);

  }
  //If your local d_min is the smallest, you've got the good itinerary, so you can print it out.
  if(d_min == rootd_min)
  {
	  // Writing a text file containing the itinerary (x,y coordinates) -
	  // to be plotted in a plotting software.
	  FILE *fp = fopen ("tsp.dat", "w");
	  for (l=0; l<N_CITIES; l++)
		{
		  fprintf (fp, "%f %f\n", x[perm_min[l]], y[perm_min[l]]);
		}
	  // Going back to the original city:
	  fprintf (fp, "%f %f\n", x[perm_min[0]], y[perm_min[0]]);
	  fclose (fp);
  }
  MPI_Finalize();


  return 0;

}
