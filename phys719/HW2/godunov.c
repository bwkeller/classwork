#include <stdio.h>
#include <math.h>
#define EPS 1e-6 //Tolerance for iterative scheme
#define GAMMA 1.4 //Adiabatic index 7/5 for diatomic fluid
#define a(s) sqrt(GAMMA*s->p/s->rho) //Macro for calculating the soundspeed a

#define TIME 1 //Total time to run for
#define NCELLS 400 //Size of the domain in cells
#define XMAX 2 //Physical size of domain
#define DX 2.0*XMAX/NCELLS //Grid spacing

//This struct defines the u state vector that the Fast Riemann Solver wants
typedef struct state
{
	float v;
	float p;
	float rho;
} state;


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

float guess_velocity(state* left, state* right)
{
	//Guess the initial flow velocity
	float z  = a(right)/a(left)*pow(left->p/right->p, (GAMMA-1)/(2*GAMMA));
	float vtl = left->v+2*a(left)/(GAMMA-1);
	float vtr = right->v-2*a(right)/(GAMMA-1);
	return (vtl*z+vtr)/(1+z);
}

void iterate(state *left, state *right, float dt)
{
	//Set up the conserved state vectors
	float rhol = left->rho;
	float ml = left->rho*left->v;
	float El = left->rho*(0.5l*left->v*left->v+left->p/(left->rho*(GAMMA-1)));
	float rhor = right->rho;
	float mr = right->rho*right->v;
	float Er = right->rho*(0.5l*right->v*right->v+right->p/(right->rho*(GAMMA-1)));

	float vstar = guess_velocity(left, right);
	float Wl, pl, plp, al; //Mach number, pressure, dp/du, sound speed for left
	float Wr, pr, prp, ar; //Mach number, pressure, dp/du, sound speed for right
	float Cl = GAMMA*left->p/a(left);
	float Cr = GAMMA*right->p/a(right);
	Wl = 0;
	Wr = 0;
	pl = 0;
	pr = 1;//set up pl and pr so that the first while comparison passes
	while(fabs(1-pl/pr) > EPS)
	{
		if (vstar <= left->v)//left-moving shock
		{
			Wl = (GAMMA+1)/4*(vstar-left->v)/a(left)-sqrt(1+pow((GAMMA+1)/4*(vstar-left->v)/a(left),2));
			pl = left->p+Cl*(vstar-left->v)*Wl;
			plp = 2*Cl*pow(Wl, 3)/(1+pow(Wl, 2));
		}
		else//left-moving rarefaction wave
		{
			al = a(left)-(GAMMA-1)/2*(vstar-left->v);
			pl = left->p*pow(al/a(left), 2*GAMMA/(GAMMA-1));
			plp = -GAMMA*pl/al;
		}
		if (vstar >= right->v)//right-moving shock
		{
			Wr = (GAMMA+1)/4*(vstar-right->v)/a(right)+sqrt(1+pow((GAMMA+1)/4*(vstar-right->v)/a(right),2));
			pr = right->p+Cr*(vstar-right->v)*Wr;
			prp = 2*Cr*pow(Wr, 3)/(1+pow(Wr, 2));
		}
		else//right-moving rarefaction wave
		{
			ar = a(right)+(GAMMA-1)/2*(vstar-right->v);
			pr = right->p*pow(ar/a(right), 2*GAMMA/(GAMMA-1));
			prp = GAMMA*pr/ar;
		}
		vstar -= (pl-pr)/(plp-prp);
	}
	float v, p, rho;
	float m, E; //specific momentum & energy
	if (vstar <= left->v)//left-moving shock
	{
		al = a(left)*sqrt(((GAMMA+1)+(GAMMA-1)*pl/left->p)/((GAMMA+1)+(GAMMA-1)*left->p/pl));
		v = left->v+a(left)*Wl;
		p = pl;
		rho =  GAMMA*pl/pow(al, 2);
		m = v*rho;
		E = 0.5*v*v+p/(rho*(GAMMA-1));
		rhol += dt/DX*m;
		ml += dt/DX*(v*m+p);
		El += dt/DX*(E+p)*v;
		rhor -= dt/DX*m;
		mr -= dt/DX*(v*m+p);
		Er -= dt/DX*(E+p)*v;
	}
	else//left-moving rarefaction wave
	{
		v = 0.5*(left->v-a(left) + vstar-al);//Average the velocities
		p = pl;
		rho =  GAMMA*pl/pow(al, 2);
		m = v*rho;
		E = 0.5*v*v+p/(rho*(GAMMA-1));
		rhol += dt/DX*m;
		ml += dt/DX*(v*m+p);
		El += dt/DX*(E+p)*v;
		rhor -= dt/DX*m;
		mr -= dt/DX*(v*m+p);
		Er -= dt/DX*(E+p)*v;
	}
	if (vstar >= right->v)//right-moving shock
	{
		ar = a(right)*sqrt(((GAMMA+1)+(GAMMA-1)*pr/right->p)/((GAMMA+1)+(GAMMA-1)*right->p/pr));
		v = right->v+a(right)*Wr;
		p = pr;
		rho =  GAMMA*pr/pow(ar, 2);
		m = v*rho;
		E = 0.5*v*v+p/(rho*(GAMMA-1));
		rhol += dt/DX*m;
		ml -= dt/DX*(v*m+p);
		El += dt/DX*(E+p)*v;
		rhor -= dt/DX*m;
		mr += dt/DX*(v*m+p);
		Er -= dt/DX*(E+p)*v;
	}
	else//right-moving rarefaction wave
	{
		v = 0.5*(right->v+a(right) + vstar+ar);//Average the velocities
		/*result.right->p = pr;*/
		/*result.right->rho =  GAMMA*pr/pow(ar, 2);*/
		/*result.centre.p = pr;*/
		/*result.centre.rho =  GAMMA*pr/pow(ar, 2);*/
	}

	//Return the conserved variables to the standard stave variables
	left->rho = rhol;
	left->v = ml/rhol;
	left->p = (El-0.5*rhol*left->v*left->v)*(GAMMA-1);
	right->rho = rhor;
	right->v = mr/rhor;
	right->p = (Er-0.5*rhor*right->v*right->v)*(GAMMA-1);
}

void evolve(state domain[])
{
	int j;
	float i;
	float dt = 1;
	state leftbound = domain[0];
	state rightbound = domain[NCELLS-1];
	state tmp;
	for(i=0; i<TIME; i+=dt)
	{
		printf("Integrating, %f complete\n", (float)i/(float)TIME);
		for(j=0; j<NCELLS; j++)
		{
			/*printf("position %d ", j);*/
			if (i == 0)//Lefthand boundary
			{
				tmp = leftbound;
				iterate(&leftbound, &domain[j], dt);
				leftbound = tmp;//Keep the boundary conditions fixed
			}
			else if(i == NCELLS-1)
			{
				tmp = rightbound;
				iterate(&domain[j], &rightbound, dt);
				rightbound = tmp;//Keep the boundary conditions fixed
			}
			else
			{
				iterate(&domain[j], &domain[j+1], dt);
			}
		}
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
