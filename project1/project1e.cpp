/*
To compile this code, run the following command:
make -f makefile_pj1e
*/

#include <armadillo>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <ctime>
using namespace arma;
using namespace std;

// Functions as used in tasks a)-d)
inline double f_mark__(double x, double h) {return pow(h,2)*100.0*exp(-10.0*x);}
inline double u__(double x)	{return (1 - (1 - exp(-10))*x - exp(-10*x));}

int main(int argc, char* argv[])
{
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
  mat A = zeros(n,n);   // LHS matrix
  double h = 1./(n+1);

// Iterate through rows filling fill main diagonal and adjacent values
  for (int i=0; i < n-1; i++)
  {
    A(i, i)    =  2;      // Main diagonal
    A(i, i+1)  = -1;      // Upper diagonal values
    A(i+1, i)  = -1;      // Lower diagonal values
  }

// Set 2nd last main diagonal manually
  A(n-1, n-1) = 2;

  mat f_mark = zeros(n);    // m x 1 matrix (vector)
  double *x = new double[n];
  double *u = new double[n];

  x[0] = h;
  f_mark(0) = f_mark__(x[0], h);
  for (int i=1; i < n ; i++)		// create x-array and compute RHS
  {
    x[i] = x[i-1] + h;
    f_mark(i) = f_mark__(x[i],h);
  }

  // Initialize variables to calculate execution time
  clock_t t;
  t = clock();

// Solve using LU decomposition
  vec v = solve(A, f_mark);

  // Calculate time by using number of clock ticks elapsed
  t = clock() - t;

  double total_ms;
  total_ms = 1000*float(t)/CLOCKS_PER_SEC;	// num. of seconds algorithm takes to run
  /* Write execution time for main algorithm to file. Used later for calculating average CPU time to compute main algorithm.*/
  ofstream ofile;
  ofile.open("CPUtime.txt", ofstream::app);
  ofile << n << " " << total_ms << endl;
  ofile.close();

  printf ("CPU time for main algorithm: %g ms\n", total_ms);

// Write exact solution and numerical solution as function of x
  ofstream myfile;
  char *project1_e_data;
  myfile.open ("project1_e_data.txt");
  for (int i=0; i < n; i++)
  {
    u[i] = u__(x[i]);						// analytical solution
    myfile << x[i] << ' ' << u[i] << ' ' << v[i] << endl;
  }

  myfile.close();
  printf ("Solution computed for n = %i, using Armadillo library. Results written to 'project1_e_data.txt'\n", n);

  return 0;
}
