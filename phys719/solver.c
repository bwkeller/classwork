#include <stdio.h>
#include <math.h>
#define EPS 1e-6 //Tolerance for iterative scheme
#define GAMMA 1.4 //Adiabatic index 7/5 for diatomic fluid
#define a(s) sqrt(GAMMA*s.p/s.rho) //Macro for calculating the soundspeed a

//This struct defines the u state vector.
typedef struct state
{
	float v;
	float p;
	float rho;
} state;

int read_input(state* left, state* right, char filename[])
{
	/* Input file format should store the ICs in the format:
	 * v_l v_r
	 * p_l p_r
	 * rho_l rho_r
	 */
	FILE * infile = fopen(filename, "r+");
	if(infile == NULL)
	{
		printf("Error reading input file");
		return 1;
	}
	fscanf(infile, "%f %f", &left->v, &right->v);
	fscanf(infile, "%f %f", &left->p, &right->p);
	fscanf(infile, "%f %f", &left->rho, &right->rho);
	fclose(infile);
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

void iterate(state left, state right)
{
	float vstar = guess_velocity(left, right);
	float Wl, pl, plp, al; //Mach number, pressure, dp/du, sound speed for left
	float Wr, pr, prp, ar; //Mach number, pressure, dp/du, sound speed for right
	float Cl = GAMMA*left.p/a(left);
	float Cr = GAMMA*right.p/a(right);
	Wl = 0;
	Wr = 0;
	al = a(left);
	ar = a(right);
	pl = 0;
	pr = 1;//set up pl and pr so that the first while comparison passes
	printf("\nInitial Guess: %5.3f\n", vstar);
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
			plp = -1*GAMMA*pl/al;
		}
		if (vstar >= right.v)//right-moving shock
		{
			Wr = (GAMMA+1)/4*(vstar-right.v)/a(right)+sqrt(1+pow((GAMMA+1)/4*(vstar-right.v)/a(right),2));
			pr = right.p+Cr*(vstar-right.v)*Wr;
			prp = 2*Cl*pow(Wr, 3)/(1+pow(Wr, 2));
		}
		else//right-moving rarefaction wave
		{
			ar = a(right)+(GAMMA-1)/2*(vstar-right.v);
			pr = right.p*pow(ar/a(right), 2*GAMMA/(GAMMA-1));
			prp = GAMMA*pr/ar;
		}
		vstar -= (pl-pr)/(plp-prp);
	}
	printf("\nFinal State\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("u*:%5.3f\n", vstar);
	printf("Side\tType\tV( head & tail)\tP\trho\n");
	if (vstar <= left.v)//left-moving shock
	{
		al = a(left)*sqrt(((GAMMA+1)+(GAMMA-1)*pl/left.p)/((GAMMA+1)+(GAMMA-1)*left.p/pl));
		printf("left\tShock\t%5.3f\t\t%5.3f\t%5.3f\n", left.v+al*Wl, pl, GAMMA*pl/pow(al, 2));
	}
	else//left-moving rarefaction wave
	{
		printf("left\tFan\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", left.v-a(left), vstar-al, pl, GAMMA*pl/pow(al, 2));
	}
	if (vstar >= right.v)//right-moving shock
	{
		ar = a(right)*sqrt(((GAMMA+1)+(GAMMA-1)*pr/right.p)/((GAMMA+1)+(GAMMA-1)*right.p/pr));
		printf("right\tShock\t%5.3f\t\t%5.3f\t%5.3f\n", right.v+ar*Wr, pr, GAMMA*pr/pow(ar, 2));
	}
	else//right-moving rarefaction wave
	{
		printf("right\tFan\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", right.v+a(right), vstar+ar, pr, GAMMA*pr/pow(ar, 2));
	}

}
int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printf("Too few arguments.  Please give me an initial condition file\n");
		return 1;
	}
	state left;
	state right;
	int success = read_input(&left, &right, argv[1]);
	if (success) return 2;
	printf("Starting Iteration: Gamma = %3.2f Eps = %3.2e\n", GAMMA, EPS);
	printf("Initial State\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("Side\tV\tP\trho\n");
	printf("left\t%5.3f\t%5.3f\t%5.3f\n", left.v, left.p, left.rho);
	printf("right\t%5.3f\t%5.3f\t%5.3f\n", right.v, right.p, right.rho);
	iterate(left, right);
	return 0;
}
