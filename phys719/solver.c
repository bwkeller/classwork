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

void read_input(state* left, state* right)
{
	/* Input file format should store the ICs in the format:
	 * v_l v_r
	 * p_l p_r
	 * rho_l rho_r
	 */
	FILE * infile = fopen("input.dat", "r");
	fscanf(infile, "%f", left->v);
	fscanf(infile, "%f", right->v);
	fscanf(infile, "%f", left->p);
	fscanf(infile, "%f", right->p);
	fscanf(infile, "%f", left->rho);
	fscanf(infile, "%f", right->rho);
	fclose(infile);
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
	pl = 0;
	pr = 1;//set up pl and pr so that the first while comparison passes
	while(abs(1-pl/pr) > EPS)
	{
		if (vstar <= left.v)
		{
			Wl = (GAMMA+1)/4*(vstar-left.v)/a(left)-sqrt(1+pow((GAMMA+1)/4*(vstar-left.v)/a(left),2));
			pl = left.p+Cl*(vstar-left.v)*Wl;
			plp = 2*Cl*pow(Wl, 3)/(1+pow(Wl, 2));
		}
		else
		{
			al = a(left)-(GAMMA-1)/2*(vstar-left.v);
			pl = left.p*pow(al/a(left), 2*GAMMA*(GAMMA-1));
			plp = -1*GAMMA*pl/al;
		}
		if (vstar >= right.v)
		{
			Wr = (GAMMA+1)/4*(vstar-right.v)/a(right)+sqrt(1+pow((GAMMA+1)/4*(vstar-right.v)/a(right),2));
			pr = right.p+Cr*(vstar-right.v)*Wl;
			prp = 2*Cl*pow(Wr, 3)/(1+pow(Wr, 2));
		}
		else
		{
			ar = a(right)-(GAMMA-1)/2*(vstar-right.v);
			pr = right.p*pow(ar/a(right), 2*GAMMA*(GAMMA-1));
			prp = -1*GAMMA*pr/ar;
		}
		vstar -= (pl-pr)/(plp-prp);
		printf("%f\t%f\t%f\n", left.v, vstar, right.v);
	}
}
int main(int argc, char* argv[])
{
	state left;
	state right;
	read_input(&left, &right);

	return 0;
}
