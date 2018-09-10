#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <iomanip>

/*
Task c) With identical values along the diagonal, we can set all elements in vectors a and c to -1.
*/

//f_mark = h**2 f
inline double f_mark__(double x, double h) {return pow(h,2)*100.0*exp(-10.0*x);}
inline double u__(double x)	{return (1 - (1 - exp(-10))*x - exp(-10*x));}
/*
Function that creates an array of size 'm' and allocates memory for indices of type double
*/
double* vec_mem(int m)
{
	double *v = new double[m];
	return v;
}
/*
Integer version of function above
*/
int* vec_mem_int(int m)
{
	int *v = new int[m];
	return v;
}

using namespace std;

int main(int argc, char* argv[])
{
	// Check if number of command line arguments are correct
	if (argc < 2 )
	{
		cout << "Input error: specify number of grid points 'n' as an additional argument." << endl;
		exit(1);
	}
	if (argc > 2)
	{
		cout << "Input error: too many arguments. Only specify number of grid points 'n'." << endl;
		exit(1);
	}
		int n=atoi(argv[1]);
		double h = 1. / (n+1);

		// Allocate memory to arrays used
		double *x, *f_mark, *u, *b, *v, *f_tilde, *b_tilde;
		x = vec_mem(n+2); f_mark = vec_mem(n+2); u = vec_mem(n+2);
		b = vec_mem(n+2); f_tilde = vec_mem(n+2); b_tilde = vec_mem(n+2);
		v = vec_mem(n+2);															// unknown solution

		for (int i=0; i < n+2 ; i++)		// create x-array
		{
			x[i] = h * i;
			f_mark[i] = f_mark__(x[i], h);
			u[i] = u__(x[i]);						// analytical solution
		}
		// Set boundary conditions
		v[0] = 0;
		v[n+1] = 0;

		// Fill main diagonal
		for (int i=0; i < n+2; i++)
		{
			b[i] = 2;
		}

		// Set initial f_tilde and b_tilde values to initial f and b values, since
		// tilde values are not used until index i = 1
		f_tilde[0] = f_mark[0];
		f_tilde[1] = f_mark[1];
		b_tilde[0] = b[0];
		b_tilde[1] = b[0];

		// Initialize variables to calculate execution time
		clock_t t;
		t = clock();

		// Main algorithm: perform forward and backward substitution
		double _b_tilde;
		for (int i=2; i < n+2; i++)
		{
			_b_tilde = 1./(b_tilde[i-1]);	// save FLOPS by calculating once
			b_tilde[i] = b[i] - _b_tilde;
	    f_tilde[i] = f_mark[i] + f_tilde[i-1]*_b_tilde;
		}
		for (int i=n; i > 0; i--)
		{
			v[i] = (f_tilde[i] + v[i+1])/b_tilde[i];
		}

		// Calculate time by using number of clock ticks elapsed
		t = clock() - t;
		double total_ms;
		total_ms = 1000*float(t)/CLOCKS_PER_SEC;	// num. of milliseconds algorithm takes to run

		/* Write execution time for main algorithm to file. Used later for calculating average CPU time to compute main algorithm.*/
		ofstream ofile;
		ofile.open("error.txt", ofstream::app);
		ofile << n << " " << total_ms << endl;
		ofile.close();
		
		printf ("CPU time for main algorithm: %g ms\n", total_ms);

		// Write exact solution and numerical solution as function of x
		ofstream myfile;
		char *project1_c_data;
		myfile.open ("project1_c_data.txt");

		for (int i=0; i < n+2; i++)
		{
			myfile << x[i] << ' ' << u[i] << ' ' << v[i] << endl;
		}
		myfile.close();
		printf ("Solution computed for n = %i. Results written to file 'project1_c_data.txt'.\n", n);

		return 0;
}
