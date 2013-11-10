#include <stdio.h>
#include <math.h>
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
	float z  = a(right)/a(left)*pow(left.p/right.p, (GAMMA-1)/(2*GAMMA));
	return z;
}

int main(int argc, char* argv[])
{
	state left;
	state right;
	read_input(&left, &right);

	return 0;
}
