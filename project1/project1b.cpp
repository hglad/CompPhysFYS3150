#define _USE_MATH_DEFINES

#include <iostream>
//#include <armadillo>
#include <cmath>
#include <fstream>

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
		int n=atoi(argv[1]);		   // endpoint of x-array
		double h = 1. / (n+1);

		// Allocate memory to arrays used
		double *x, *f_mark, *u, *b, *v, *f_tilde, *b_tilde;
		int *a, *c;
		x = vec_mem(n+2); f_mark = vec_mem(n+2); u = vec_mem(n+2);
		b = vec_mem(n+2); f_tilde = vec_mem(n+2); b_tilde = vec_mem(n+2);
		a = vec_mem_int(n+1); c = vec_mem_int(n+1);		// lower and upper diagonals
		v = vec_mem(n+2);															// unknown solution

		for (int i=0; i < n+2 ; i++)		// create x-array
		{
			x[i] = h * i;
			f_mark[i] = f_mark__(x[i], h);
			u[i] = u__(x[i]);						// analytical solution
		//	printf ("x = %g, f_mark = %g \n", x[i], f_mark[i]);
		}
		// Set boundary conditions
		v[0] = 0;
		v[n+1] = 0;

		// Fill lower and upper diagonals
		for (int i=0; i < n+1; i++)
		{
			a[i] = -1;
			c[i] = -1;
		}
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
		double a_b_tilde;
		for (int i=2; i < n+2; i++)
		{
			a_b_tilde = a[i-1]/b_tilde[i-1];	// save FLOPS by only calculating once
			b_tilde[i] = b[i] - a_b_tilde * c[i-1];
	    f_tilde[i] = f_mark[i] - a_b_tilde * f_tilde[i-1];
		}
		for (int i=n; i > 0; i--)
		{
			v[i] = (f_tilde[i] - v[i+1]*c[i])/b_tilde[i];
		}

		// Calculate time by using number of clock ticks elapsed
		t = clock() - t;
		double total_seconds;
		total_seconds = float(t)/CLOCKS_PER_SEC;	// num. of seconds algorithm takes to run
		printf ("CPU time for main algorithm: %g seconds\n", total_seconds);

		// Write arrays x, u and v to file
		ofstream myfile;
		char *project1_b_data;
		myfile.open ("project1_b_data.txt");

		for (int i=0; i < n+2; i++)
		{
		//	cout << i << ' ' << b_tilde[i] << endl;
		//	printf ("%g %g %g\n", x[i], f_tilde[i], v[i]);
			myfile << x[i] << ' ' << u[i] << ' ' << v[i] << endl;
		}
		cout << v[n+1] << endl;
		myfile.close();
		printf ("Solution computed for n = %i. Results written to file.\n", n);

		return 0;
}
