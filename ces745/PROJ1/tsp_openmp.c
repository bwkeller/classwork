/* Exercise to solve the Traveling Salesman Problem (TSP) using
 a brute force Monte Carlo approach. Your task is to convert the serial
 computation to a parallel computation. The code has to print the
 current shortest total traveling distance found fairly frequently, so 
 you can monitor its progress.

 We use random number generator rand_r() specifically because it is thread-safe
 (it doesn't keep anything in a global state; everything is defined by the 
 variable "seed", which can be made private to the thread). 

 You only have to parallelize the code located between the
 "----------" comment lines. So e.g. don't bother to compute the
 distance matrix in parallel - keep it as it is. The reason for
 that is the distance matrix calculations scale as N^2, whereas the
 minimum traveling distance calculations scale as (N-1)! (N is the number
 of cities), so for N>6 the former becomes totally negligible compared
 to the latter in terms of CPU cycles.

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

 - Try to minimize the impact of a critical region as much as possible
 (OpenMP).

 - In OpenMP, experiment with different loop schedules - you might see
 big difference. Find the most efficient kind of schedule for this setup.


To compile:

 - serial:
icc -O2 tsp.c -o tsp

 - OpenMP:
icc -O2 -openmp tsp_omp.c -o tsp_omp


*/

#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>

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
  struct timeval  tdr0, tdr1, tdr;
  double restime;
  int devid, devcount, error, i, j, k, l, l1, ltemp;
  int perm[N_CITIES], perm_min[N_CITIES];
  float d;

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

  float d_min = 1e30;
  unsigned int seed = 111;
#pragma omp parallel default(none) shared(dist, d_min, perm_min) private(perm, k, l, d, l1, ltemp, restime, tdr1, tdr) firstprivate(seed, tdr0)
  {
	  // We need each thread to have a unique RNG, so that each thread actually tests different paths
	  // Adding the thread number to the initial seed is an easy way to do this.
	  seed += omp_get_thread_num();
	  // We always start from the same (0-th) city:
	  perm[0] = 0;
	  

	  // Cycle for Monte Carlo steps:
	  //
	  // Use the guided scheduler here, it is *slightly* faster than the static
	  // scheduler, and almost 10x faster than dynamic, which is horrendously slow
#pragma omp for schedule(guided)
	  for (k=0; k<N_MC; k++)
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

		  // Finding globally shortest TSP distance:
		  if (d < d_min)
		{
#pragma omp critical
			{
			  // Check to ensure that d_min hasn't changed since the comparison before
			  // the critical section.  Putting the critical section inside the 
			  // if(d < d_min) comparison GREATLY reduces the number of times critical
			  // is hit.
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
		}

		}  
#pragma omp single	  
	  {
		  // We don't want this timing data to be printed repeatedly by each thread
		  // so we put this section into a serial workshare loop.
		  gettimeofday (&tdr1, NULL);
		  tdr = tdr0;
		  timeval_subtract (&restime, &tdr1, &tdr);
		  printf ("Shortest total distance: %f\n", d_min);

		  printf ("%e\n", restime);
	  }
  }


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


  return 0;

}
