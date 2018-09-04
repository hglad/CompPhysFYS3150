#define _USE_MATH_DEFINES

#include <iostream>
//#include <armadillo>
#include <cmath>
#include <fstream>

//f_mark = h**2 f
inline double f_mark__(double x, double h) {return pow(h,2)*100.0*exp(-10.0*x);}
inline double u__(double x)	{return (1 - (1 - exp(-10))*x - exp(-10*x));}
using namespace std;

int main(int argc, char* argv[])
{
		int n=atoi(argv[1]);
		double x_end=1;		   // endpoint of x-array
		double x[n];			   // vector to hold x-values
		double f_mark[n]; double u[n];
		//double *p;				 // pointer to  double float
		//p = &x[0];
		double h = x_end / (n-1);

		for (int i=0; i < n ; i++)		// create x-array
		{
			x[i] = x_end * h * i;
			f_mark[i] = f_mark__(x[i], h);
			u[i] = u__(x[i]);						// analytical solution
			printf ("x = %g, f_mark = %g \n", x[i], f_mark[i]);
		}

		int a[n-1]; int c[n-1];	        // lower and upper diagonals
		double b[n];	 double v[n];							// diagonal and LHS vector
		v[0] = 0;
		v[n] = 0;

		// Fill lower and upper diagonals
		for (int i=0; i < n-1; i++)
		{
			a[i] = -1;
			c[i] = -1;
		}

		// Initialize b_tilde and f_tilde, fill with b and f_mark values initially
		double b_tilde[n]; double f_tilde[n];
		for (int i=0; i < n; i++)
		{
			b[i] = 2;
			f_tilde[i] = f_mark[i];
			b_tilde[i] = b[i];
		}
		// Main algorithm: perform forward and backward substitution
		for (int i=1; i < n; i++)
		{
			b_tilde[i] = b[i] - a[i-1]/b_tilde[i-1] * c[i-1];
	    f_tilde[i] = f_mark[i] - a[i-1]/b_tilde[i-1] * f_tilde[i-1];
			printf ("b_tilde[%i] = %g, f_tilde[%i] = %g \n", i, b_tilde[i], i, f_tilde[i]);
		}

		for (int i=n-2; i > 0; i--)
		{
			v[i] = (f_tilde[i] - v[i+1]*c[i])/b_tilde[i];
			cout << v[i] << endl;
		}

		// Write to file
		ofstream myfile;
		char *project1_b_data;
		myfile.open ("project1_b_data.txt");

		for (int i=0; i < n; i++)
		{
			myfile << x[i] << ' ' << u[i] << ' ' << v[i] << '\n' << endl;
		}
		myfile.close();

		return 0;
}
