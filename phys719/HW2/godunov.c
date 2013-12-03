#include <stdio.h>
#include <assert.h>
#include <math.h>
#define EPS 1e-6 //Tolerance for iterative scheme
#define GAMMA 1.4 //Adiabatic index 7/5 for diatomic fluid
#define a(s) sqrt(GAMMA*s.p/s.rho) //Macro for calculating the soundspeed a

#define TIME  0.01//Total time to run for
#define NCELLS 40 //Size of the domain in cells
#define XMAX 2 //Physical size of domain
#define DX (2.0*XMAX/NCELLS) //Grid spacing
#define DT (DX/4.0) //Timestep size

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

typedef struct riemann
{
    float vstar;
    int fanL; //0 for shock, 1 for fan
    int fanR; //0 for shock, 1 for fan
    state left;
    float tailL;//store the extra tail velocity for fans
    state right;
    float tailR;//store the extra tail velocity for fans
    state left0;
    state right0;
} riemann;

//Get a conservative flux from the FRS state variables
flux fluxify(state in)
{
    flux out;
    out.rhodot = in.v*in.rho;
    out.momentumdot = in.rho*in.v*in.v+in.p;
    out.energydot = in.v*(0.5*in.rho*in.v*in.v+in.p*(1+1.0/(GAMMA-1)));
    return out;
}

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


riemann iterate(state left, state right)
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
        out.left = lout;
        out.fanL = 0;
	}
	else//left-moving rarefaction wave
	{
        /*printf("%f %f\n", vstar, al);*/
        state lout = {left.v-a(left), pl, GAMMA*pl/pow(al, 2)};
        out.left = lout;
        out.tailL = vstar-al;
        out.fanL = 1;
	}
	if (vstar >= right.v)//right-moving shock
	{
		ar = a(right)*sqrt(((GAMMA+1)+(GAMMA-1)*pr/right.p)/((GAMMA+1)+(GAMMA-1)*right.p/pr));
        state rout = {right.v+a(right)*Wr, pr, GAMMA*pr/pow(ar, 2)};
        out.right = rout;
        out.fanR = 0;
	}
	else//right-moving rarefaction wave
	{
        state rout = {right.v+a(right), pr, GAMMA*pl/pow(ar, 2)};
        out.right = rout;
        out.tailR = vstar+ar;
        out.fanR = 1;
	}
    return out;

}
void update(state domain[], flux fluxes[])
{
	int i;
	for(i=1; i<NCELLS-1; i++)
	{
		//convert state variables into their conserved versions
		float rho = domain[i].rho;
		float momentum = domain[i].rho*domain[i].v;
		float energy = 0.5*domain[i].rho*domain[i].v*domain[i].v+domain[i].p/(GAMMA-1);
		//apply fluxes
		rho += DT/DX*fluxes[i].rhodot;
		momentum += DT/DX*fluxes[i].momentumdot;
		energy += DT/DX*fluxes[i].energydot;
        assert(energy > 0);
        assert(rho > 0);
		//revert to non-conserved quantities
		domain[i].v = momentum/rho;
		domain[i].rho = rho;
		domain[i].p = (energy-0.5*rho*momentum*domain[i].v)*(GAMMA-1);
        assert(domain[i].p > 0);
	}
}

flux add_flux(flux left, flux right)
{
    flux out;
    out.rhodot = left.rhodot-right.rhodot;
    out.momentumdot = left.momentumdot-right.momentumdot;
    out.energydot = left.energydot-right.energydot;
    return out;
}

flux flux_calc(riemann sol)
{
    state u0;
    if(sol.vstar < 0)//Left-moving contact
    {
        //Left moving shock or fan head
        if( (!sol.fanR && sol.right.v < 0) || (sol.fanR && sol.right.v < 0) )
        {
            u0 = sol.right0;
        }
        else if( !sol.fanR && sol.right.v > 0 )//Right-moving shock
        {
            u0 = sol.right;
        }
        else if( sol.fanR && sol.tailR > 0)//Right-moving fan tail
        {
            float vpost;//the post-fan velocity
            if(sol.fanL)
            {
                vpost = sol.vstar;
            }
            else
            {
                vpost = sol.left.v;
            }
            u0 = sol.right;
            u0.v = vpost;
        }
        else if( sol.fanR && sol.tailR < 0 && sol.right.v > 0)//Inside the fan!
        {
            float vpost;//the post-fan velocity
            if(sol.fanL)
            {
                vpost = sol.vstar;
            }
            else
            {
                vpost = sol.left.v;
            }
            float vdot = (vpost-sol.right0.v)/((sol.tailR-sol.right.v));
            float rhodot = (sol.right.rho-sol.right0.rho)/((sol.tailR-sol.right.v));
            float pdot = (sol.right.p-sol.right0.p)/((sol.tailR-sol.right.v));
            u0.v = vpost-vdot*sol.tailR;
            u0.rho = sol.right.rho-rhodot*sol.tailR;
            u0.p= sol.right.p-pdot*sol.tailR;
        }
        else
        {
            assert(0);//We should NEVER be here
        }
    }
    else//Right-moving contact
    {
        //Right moving shock or fan head
        if( (!sol.fanL && sol.left.v > 0) || (sol.fanL && sol.left.v > 0) )
        {
            u0 = sol.left0;
        }
        else if( !sol.fanL && sol.left.v < 0 )//Right-moving shock
        {
            u0 = sol.left;
        }
        else if( sol.fanL && sol.tailL < 0)//Right-moving fan tail
        {
            float vpost;//the post-fan velocity
            if(sol.fanR)
            {
                vpost = sol.vstar;
            }
            else
            {
                vpost = sol.right.v;
            }
            u0 = sol.left;
            u0.v = vpost;
        }
        else if( sol.fanL && sol.tailL > 0 && sol.left.v < 0)//Inside the fan!
        {
            float vpost;//the post-fan velocity
            if(sol.fanR)
            {
                vpost = sol.vstar;
            }
            else
            {
                vpost = sol.right.v;
            }
            float vdot = (vpost-sol.left0.v)/((sol.tailL-sol.left.v));
            float rhodot = (sol.left.rho-sol.left0.rho)/((sol.tailL-sol.left.v));
            float pdot = (sol.left.p-sol.left0.p)/((sol.tailL-sol.left.v));
            u0.v = vpost-vdot*sol.tailL;
            u0.rho = sol.left.rho-rhodot*sol.tailL;
            u0.p= sol.left.p-pdot*sol.tailL;
        }
        else
        {
            assert(0);//We should NEVER be here
        }

    }
    return fluxify(u0);
}
void evolve(state domain[])
{
	int j;
	float i;
	state tmp;
	riemann left, right;
	flux fluxes[NCELLS];
	for(i=0; i<TIME; i+=DT)
	{
		printf("Integrating, %f complete\n", (float)i/(float)TIME);
		for(j=0; j<NCELLS; j++)
		{
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
            if(j > 17 && j < 22)
            {
                flux lflux = flux_calc(left);
                flux rflux = flux_calc(right);
                /*printf("j: %d l_drho: %f l_dmom: %f l_de: %f\n", j, lflux.rhodot, lflux.momentumdot, lflux.energydot);*/
                /*printf("j: %d r_drho: %f r_dmom: %f r_de: %f\n", j, rflux.rhodot, rflux.momentumdot, rflux.energydot);*/
            }
            fluxes[j] = add_flux(flux_calc(left), flux_calc(right));
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
