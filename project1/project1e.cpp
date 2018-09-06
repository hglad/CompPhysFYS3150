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
  int n=atoi(argv[1]);
  double h = 1./(n+1);
  mat A = zeros(n,n);   // LHS matrix

/* Iterate through rows and set all diagonal elements up to n-3. Last index i = n is not iterated over since i+1 will go out of bounds at this point. */
  for (int i=0; i < n-1; i++)
  {
    A(i, i)    =  2;      // Main diagonal
    A(i, i+1)  = -1;      // Upper diagonal values
    A(i+1, i)  = -1;      // Lower diagonal values
  }
// Set last diagonal elements
  A(n-1, n-1) = 2;

  mat f_mark = zeros(n);
  double *x = new double[n];
  double *u = new double[n];

// Compute values at endpoints
  x[0] = h;
  f_mark(0) = f_mark__(x[0], h);
  x[n-1] = x[0]+(n-1)*h;
  f_mark(n-1) = f_mark__(x[n-1], h);

  for (int i=1; i < n-1 ; i++)		// create x-array and compute RHS
  {
    x[i] = x[i-1] + h;
    f_mark(i) = f_mark__(x[i],h);
  }
// Solve using LU decomposition
  vec v = solve(A, f_mark);

// Write results to file
  ofstream myfile;
  char *project1_e_data;
  myfile.open ("project1_e_data.txt");
  for (int i=0; i < n; i++)
  {
    u[i] = u__(x[i]);						// analytical solution
    myfile << x[i] << ' ' << u[i] << ' ' << v[i] << endl;
  }

  myfile.close();
  printf ("Solution computed for n = %i, using Armadillo library. Results written to file.\n", n);

  return 0;
}
