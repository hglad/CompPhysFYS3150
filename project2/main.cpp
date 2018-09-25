// Compile: g++ project2b.cpp main.cpp -o project2b.x -DARMA_DONT_USE_WRAPPER -lblas -llapack

#include "project2b.h"
using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{
  // Initialize matrix and necessary variables
  int n = atoi(argv[1]);
  double eps = 1e-8;    // tolerance to represent values close enough to 0

  mat A = zeros(n,n);
  double h = 1./(n);
  double d = 2./(h*h);
  double a = -1./(h*h);

  // Generate tridiagonal matrix
  for (int i=0; i < n-1; i++)
  {
    A(i, i)    = d;      // Main diagonal
    A(i, i+1)  = a;      // Upper diagonal values
    A(i+1, i)  = a;      // Lower diagonal values
  }
  A(n-1, n-1) = d;       // Last diagonal element

//  arma_diag_test(A);    // Find eigenvalues- and vectors with armadillo

  // Rotate matrix until value of highest element is close to 0
  cout << A << endl;       // print initial matrix

  double max_elem = 2*eps; // initialize max variable with arbitrary value > eps
  int num_rotations = 0;

  while (max_elem > eps)
  {
    Jacobi_rot(A, max_elem);
    num_rotations += 1;
    //printf("max_elem outside function: %g\n", max_elem);
  }

  cout << A << endl;    // print resulting matrix
  printf("Number of rotations performed: %i\n", num_rotations );
  return 0;

}
