#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

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
		cout << "Input error: specify which power of 10 to iterate up to." << endl;
		exit(1);
	}
	if (argc > 2)
	{
		cout << "Input error: too many arguments. Only specify power of 10 to iterate up to." << endl;
		exit(1);
	}
		ofstream myfile;
		char *project1_d_data;
		myfile.open ("project1_d_data.txt");

		int ind_counter = 0;			// counter for how many powers have been run over
		int power = atoi(argv[1]);		// Specify how many powers of 10 to run over
		for (int n=10; n <= pow(10,power); n*=10 )// Run over powers of 10
		{
			double h = 1. / (n+1);

			// Allocate memory to arrays used
			double *x, *f_mark, *u, *b, *v, *f_tilde, *b_tilde;
			double *eps;			// error estimate array
			x = vec_mem(n+2); f_mark = vec_mem(n+2); u = vec_mem(n+2);
			b = vec_mem(n+2); f_tilde = vec_mem(n+2); b_tilde = vec_mem(n+2);
			v = vec_mem(n+2);															// unknown solution
			eps = vec_mem(n+2);

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

			/* Find max error by running over all points and checking if current error is larger than previous error, then add to list containing max relative error values */
			double max_error[power];
			for (int i=1; i < n+1; i++)
			{
				eps[i] = log10(abs( (v[i] - u[i])/u[i] ));
				if ( eps[i] > eps[i-1])
				{
					max_error[ind_counter] = eps[i];
				}
			}

			printf ("Solution computed for h = %g. Max relative error: %g\n", h, max_error[ind_counter]);

			// Write maximum error for current step size to file
			myfile << log10(h) << ' ' << max_error[ind_counter] << endl;
			ind_counter += 1;
		}

		myfile.close();
		printf ("Results written to file 'project1_d_data.txt'.\n");
		return 0;
}
