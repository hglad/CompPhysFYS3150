#include "project2.h"

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]);
  double rho_max = atoi(argv[2]);
  double rho_0 = 0;
  double h = rho_max/(n);    // rho_min = 0
  double d = 2./(h*h);
  double a = -1./(h*h);

  vec rho = linspace(rho_0 + h, rho_max - h, n);
  mat A = zeros(n,n);
  // Generate tridiagonal matrix with added potential
  for (int i=0; i < n-1; i++)
  {
    A(i, i)    = d + rho[i]*rho[i]; // Main diagonal with potential
    A(i, i+1)  = a;                 // Upper diagonal values
    A(i+1, i)  = a;                 // Lower diagonal values
  }
  A(n-1, n-1) = d + rho[n-1]*rho[n-1];       // Last diagonal element

  //cout << A << endl;
  double eps = 1e-8;
  double max_elem = 2*eps; // initialize max variable with arbitrary value > eps
  int num_rotations = 0;

  while (max_elem > eps)
  {
    Jacobi_rot(A, max_elem);
    num_rotations += 1;
    //printf("max_elem outside function: %g\n", max_elem);
  }
  vec eigvals = A.diag();
  cout << eigvals << endl;

  printf("Number of rotations performed: %i\n", num_rotations );

  return 0;
}
