// Compile: g++ project2.cpp main_b.cpp -o project2b.x -DARMA_DONT_USE_WRAPPER -lblas -llapack
#include "project2.h"

int main(int argc, char *argv[])
{
  // Initialize matrix and necessary variables
  int n = atoi(argv[1]);

  double h = 1./(n);
  vec d = 2./(h*h) * ones<vec>(n);
  double a = -1./(h*h);

  // Generate tridiagonal matrix
  mat A = generate_mat(n, d, a, a);
  vec eigvals;
  eigvals = arma_diag(A);    // Find eigenvalues- and vectors with armadillo
  cout << eigvals << endl;
  // Rotate matrix until value of highest element is close to 0
  cout << A << endl;       // print initial matrix

  double eps = 1e-8;    // tolerance to represent values close enough to 0
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
