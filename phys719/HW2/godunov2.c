#include <stdio.h>
#include <assert.h>
#include <math.h>
#define EPS 1e-6 //Tolerance for iterative scheme
#define GAMMA 1.4 //Adiabatic index 7/5 for diatomic fluid
#define a(s) sqrt(GAMMA*s.p/s.rho) //Macro for calculating the soundspeed a

#define TIME  1e-5//Total time to run for
#define NCELLS 4 //Size of the domain in cells
#define XMAX 2 //Physical size of domain
#define DX (2.0*XMAX/NCELLS) //Grid spacing
#define DT (DX/200.0) //Timestep size

//This struct defines the u state vector that the Fast Riemann Solver wants
typedef struct state
{
	float v;
	float p;
	float rho;
} state;

typedef struct riemann
{
    float vstar;
    state left;
    float tailL;//store the extra tail velocity for fans
    state right;
    float tailR;//store the extra tail velocity for fans
    state left0;
    state right0;
} riemann;

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
		domain[i].v += 0.08;
		/*domain[i].v *= -1;*/
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


riemann solve_riemann(state left, state right)
{
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
    riemann out;
    out.vstar = vstar;
    out.left0 = left;
    out.right0 = right;
	if (vstar <= left.v)//left-moving shock
	{
		al = a(left)*sqrt(((GAMMA+1)+(GAMMA-1)*pl/left.p)/((GAMMA+1)+(GAMMA-1)*left.p/pl));
        state lout = {left.v+a(left)*Wl, pl, GAMMA*pl/pow(al, 2)};
        out.tailL = left.v+a(left)*Wl;
        out.left = lout;
	}
	else//left-moving rarefaction wave
	{
        state lout = {left.v-a(left), pl, GAMMA*pl/pow(al, 2)};
        out.left = lout;
        out.tailL = vstar-al;
	}
	if (vstar >= right.v)//right-moving shock
	{
		ar = a(right)*sqrt(((GAMMA+1)+(GAMMA-1)*pr/right.p)/((GAMMA+1)+(GAMMA-1)*right.p/pr));
        state rout = {right.v+a(right)*Wr, pr, GAMMA*pr/pow(ar, 2)};
        out.tailR = right.v+a(right)*Wr;
        out.right = rout;
	}
	else//right-moving rarefaction wave
	{
        state rout = {right.v+a(right), pr, GAMMA*pl/pow(ar, 2)};
        out.right = rout;
        out.tailR = vstar+ar;
	}
    return out;

}
state interpolate(riemann interface, int side)
{
	state out;
	float vmax, vmin, rhomax, rhomin, pmax, pmin;
	if(side)//use left wave
	{
		vmin = interface.left0.v;
		pmin = interface.left0.p;
		rhomin = interface.left0.rho;
		pmax = interface.left.p;
		rhomax = interface.left.rho;
		if(interface.tailR == interface.right.v)//right wave's a shock
		{
			vmax = interface.right.v;
		}
		else//right wave's a rarefaction
		{
			vmax = interface.vstar;
		}
		out.v = vmin-(vmax-vmin)/(interface.tailL-interface.left.v)*interface.left.v;
		out.p = pmin-(pmax-pmin)/(interface.tailL-interface.left.v)*interface.left.v;
		out.rho = rhomin-(rhomax-rhomin)/(interface.tailL-interface.left.v)*interface.left.v;
	}
	else
	{
		vmin = interface.right0.v;
		pmin = interface.right0.p;
		rhomin = interface.right0.rho;
		pmax = interface.right.p;
		rhomax = interface.right.rho;
		if(interface.tailL == interface.left.v)//left wave's a shock
		{
			vmax = interface.left.v;
		}
		else//left wave's a rarefaction
		{
			vmax = interface.vstar;
		}
		out.v = vmin-(vmax-vmin)/(interface.tailL-interface.right.v)*interface.right.v;
		out.p = pmin-(pmax-pmin)/(interface.tailL-interface.right.v)*interface.right.v;
		out.rho = rhomin-(rhomax-rhomin)/(interface.tailL-interface.right.v)*interface.right.v;
	}
	return out;
}
state get_interface_state(riemann interface)
{
	if(interface.right0.v == interface.left0.v && 
			interface.right0.p == interface.left0.p &&
			interface.right0.rho == interface.right0.rho)//weird equal-states case
	{
		printf("Equal!\n");
		return interface.right0;
	}
	if(interface.right.v < 0)//Right wave has moved past boundary
	{
		printf("6\n");
		return interface.right0;
	}
	if(interface.tailR < 0 && interface.tailR != interface.right.v)//Inside a right-moving fan
	{
		printf("5\n");
		return interpolate(interface, 0);//get the right-fan solution
	}
	if(interface.vstar < 0)//Between the right-wave and contact
	{
		printf("4");
		if(interface.tailR == interface.right.v)//right wave is a shock
		{
			printf(" shock\n");
			return interface.right;
		}
		else//right wave is a fan
		{
			state u0 = interface.right;
			if(interface.tailR == interface.right.v)//left wave is a fan
			{
				u0.v = interface.vstar;//Use the contact speed
			}
			else
			{
				u0.v = interface.left.v;//Use the left shock speed
			}
			printf(" fan\n");
			return u0;
		}
	}
	if(interface.tailL < 0)//Between the contact and the left-wave
	{
		printf("3");
		if(interface.tailL == interface.left.v)//left wave is a shock
		{
			printf(" shock\n");
			return interface.left;
		}
		else//left wave is a fan
		{
			state u0 = interface.left;
			if(interface.tailL == interface.left.v)//right wave is a fan
			{
				u0.v = interface.vstar;//Use the contact speed
			}
			else
			{
				u0.v = interface.right.v;//Use the right shock speed
			}
			printf(" fan\n");
			return u0;
		}
	}
	if(interface.left.v < 0 && interface.tailL != interface.left.v)//Inside a left-moving fan
	{
		printf("2\n");
		return interpolate(interface, 1);//get the left-fan solution
	}
	else//Left wave has moved past boundary
	{
		printf("1\n");
		return interface.left0;
	}
}
//Calculate & Apply fluxes over the range t={0,TIME}
void evolve(state domain[])
{
	int j;
	float i;
	state u0;//the state at the interface
	riemann interface;
	for(i=0; i<TIME; i+=DT)
	{
		printf("Integrating, %f complete\n", i/(float)TIME);
		for(j=0; j<NCELLS-1; j++)
		{
			interface = solve_riemann(domain[j], domain[j+1]);
			u0 = get_interface_state(interface);
			printf("%d %f %f %f\n", j, u0.v, u0.p, u0.rho);
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
