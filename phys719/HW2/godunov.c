#include <stdio.h>
#include <math.h>
#define EPS 1e-6 //Tolerance for iterative scheme
#define GAMMA 1.4 //Adiabatic index 7/5 for diatomic fluid
#define a(s) sqrt(GAMMA*s.p/s.rho) //Macro for calculating the soundspeed a

#define TIME 1 //Total time to run for
#define NCELLS 400 //Size of the domain in cells
#define XMAX 2 //Physical size of domain
#define DX 2.0*XMAX/NCELLS //Grid spacing
/*#define DT DX/2 //Timestep size*/
#define DT 2

//This struct defines the u state vector that the Fast Riemann Solver wants
typedef struct state
{
	float v;
	float p;
	float rho;
} state;

//The conserved variables (density, momentum, energy)
typedef struct flux
{
	float rhodot;
	float momentumdot;
	float energydot;
} flux;

int read_input(state domain[], char filename[])
{
	/* Input file format should store the ICs in the format:
	 * v p rho
	 * for each cell in the grid
	 */
	FILE * infile = fopen(filename, "r+");
	if(infile == NULL)
	{
		printf("Error reading input file");
		return 1;
	}
	int i;
	for(i=0; i<NCELLS; i++)
	{
		fscanf(infile, "%f %f %f", &(domain[i].v), &(domain[i].p), &(domain[i].rho));
	}
	fclose(infile);
	return 0;
}

int write_output(state domain[])
{
	FILE * outfile = fopen("output.dat", "w+");
	if(outfile == NULL)
	{
		printf("Error opening output file");
		return 1;
	}
	int i;
	for(i=0; i<NCELLS; i++)
	{
		fprintf(outfile, "%f %f %f\n", domain[i].v, domain[i].p, domain[i].rho);
	}
	fclose(outfile);
	return 0;
}

float guess_velocity(state left, state right)
{
	//Guess the initial flow velocity
	float z  = a(right)/a(left)*pow(left.p/right.p, (GAMMA-1)/(2*GAMMA));
	float vtl = left.v+2*a(left)/(GAMMA-1);
	float vtr = right.v-2*a(right)/(GAMMA-1);
	return (vtl*z+vtr)/(1+z);
}

flux iterate(state left, state right)
{
	flux output = {0,0,0};
	float vstar = guess_velocity(left, right);
	float Wl, pl, plp, al; //Mach number, pressure, dp/du, sound speed for left
	float Wr, pr, prp, ar; //Mach number, pressure, dp/du, sound speed for right
	float Cl = GAMMA*left.p/a(left);
	float Cr = GAMMA*right.p/a(right);
	Wl = 0;
	Wr = 0;
	pl = 0;
	pr = 1;//set up pl and pr so that the first while comparison passes
	while(fabs(1-pl/pr) > EPS)
	{
		if (vstar <= left.v)//left-moving shock
		{
			Wl = (GAMMA+1)/4*(vstar-left.v)/a(left)-sqrt(1+pow((GAMMA+1)/4*(vstar-left.v)/a(left),2));
			pl = left.p+Cl*(vstar-left.v)*Wl;
			plp = 2*Cl*pow(Wl, 3)/(1+pow(Wl, 2));
		}
		else//left-moving rarefaction wave
		{
			al = a(left)-(GAMMA-1)/2*(vstar-left.v);
			pl = left.p*pow(al/a(left), 2*GAMMA/(GAMMA-1));
			plp = -GAMMA*pl/al;
		}
		if (vstar >= right.v)//right-moving shock
		{
			Wr = (GAMMA+1)/4*(vstar-right.v)/a(right)+sqrt(1+pow((GAMMA+1)/4*(vstar-right.v)/a(right),2));
			pr = right.p+Cr*(vstar-right.v)*Wr;
			prp = 2*Cr*pow(Wr, 3)/(1+pow(Wr, 2));
		}
		else//right-moving rarefaction wave
		{
			ar = a(right)+(GAMMA-1)/2*(vstar-right.v);
			pr = right.p*pow(ar/a(right), 2*GAMMA/(GAMMA-1));
			prp = GAMMA*pr/ar;
		}
		vstar -= (pl-pr)/(plp-prp);
	}

	return output;
}

void update(state domain[], flux fluxes[])
{
	int i;
	for(i=0; i<NCELLS; i++)
	{
		//convert state variables into their conserved versions
		float rho = domain[i].rho;
		float momentum = domain[i].rho*domain[i].v;
		float energy = 0.5*domain[i].rho*domain[i].v*domain[i].v+domain[i].p/(GAMMA-1);
		//apply fluxes
		rho += DT/DX*fluxes[i].rhodot;
		momentum += DT/DX*fluxes[i].momentumdot;
		energy += DT/DX*fluxes[i].energydot;
		//revert to non-conserved quantities
		domain[i].v = momentum/rho;
		domain[i].rho = rho;
		domain[i].p = (energy-0.5*rho*momentum*domain[i].v)*(GAMMA-1);
	}
}

void evolve(state domain[])
{
	int j;
	float i;
	state tmp;
	flux left, right;
	flux fluxes[NCELLS];
	for(i=0; i<TIME; i+=DT)
	{
		printf("Integrating, %f complete\n", (float)i/(float)TIME);
		for(j=0; j<NCELLS; j++)
		{
			/*printf("position %d ", j);*/
			if (j == 0)//Lefthand boundary
			{
				left = iterate(domain[j], domain[j]);
				right = iterate(domain[j], domain[j+1]);
			}
			else if(j == NCELLS-1)//Righthand boundary
			{
				left = iterate(domain[j-1], domain[j]);
				right = iterate(domain[j], domain[j]);
			}
			else
			{
				left = iterate(domain[j-1], domain[j]);
				right = iterate(domain[j], domain[j+1]);
			}
			fluxes[j].rhodot = left.rhodot-right.rhodot;
			fluxes[j].momentumdot = left.momentumdot-right.momentumdot;
			fluxes[j].energydot = left.energydot-right.energydot;
		}
		update(domain, fluxes);
	}
}

int main(int argc, char* argv[])
{
	state domain[NCELLS];
	if(argc < 2)
	{
		printf("Too few arguments.  Please give me an initial condition file\n");
		return 1;
	}
	int success = read_input(domain, argv[1]);
	if (success) return 2;
	evolve(domain);
	success = write_output(domain);
	if (success) return 2;
	return 0;
}
