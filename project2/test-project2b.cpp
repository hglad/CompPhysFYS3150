// Compile: g++ project2b.cpp test-project2b.cpp -o test.x -DARMA_DONT_USE_WRAPPER -lblas -llapack

// Test that our method finds the eigenvalues of a random 5x5 matrix
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "project2b.h"

vec jacobi_diag(mat& B)
{
  int m = B.n_cols;     // dimension of matrix
  double eps = 1e-8;
  double max_elem = 2*eps; // initialize max variable with arbitrary value > eps
  int num_rotations = 0;
  vec eigvals;

  while (max_elem > eps)
  {
    Jacobi_rot(B, max_elem);
    num_rotations += 1;
  }
  // Store eigenvalues in vector
  eigvals = B.diag();

  cout << sort(eigvals) << endl;

  // Return eigenvalues sorted from lowest to highest in value (like Armadillo does)
  return sort(eigvals);
}

vec arma_diag(mat B)
{
  // Generate vector to hold eigenvalues and matrix to hold eigenvectors
  vec eigval;
  mat eigvec;
  // Decompose into eigenvalues and eigenvalues
  eig_sym(eigval, eigvec, B);
  cout << "Eigenvalues and eigenvectors:" << endl;
  cout << eigval << eigvec << endl;
  return eigval;
}

mat generate_mat(int m)
{
  mat B = zeros(m,m);

  srand (time(NULL));
  double d = rand() % 1000;
  double a = rand() % 1000;
  // Test for a tridiagonal matrix. Dense matrix produces less accurate results.

  for (int i=0; i < m-1; i++)
  {
    B(i, i)    = d;      // Main diagonal
    B(i, i+1)  = a;      // Upper diagonal values
    B(i+1, i)  = -a;      // Lower diagonal values
  }

  B(m-1, m-1) = d;       // Last diagonal element
  return B;
}

int m = 5;
mat B = generate_mat(m);

TEST_CASE("Check if eigenvalues are correct")
{
  vec eig_jac = jacobi_diag(B);  // Find eigenvalues using Jacobi's method
  vec eig_arma = arma_diag(B);   // Generate eigenvalue vector with Armadillo
  for (int i=0; i < m; i++)
  {
    REQUIRE( eig_jac(i) == Approx(eig_arma(i)) );
  }
}

TEST_CASE("Check if our method finds the maximum offdiagonal element")
{
  B.diag().fill(0);
  double max_arma = B.max();
  double max_elem;
  int k, l;
  max_elem = max_offdiag(B, &k, &l, max_elem);

  cout << max_arma << ' ' << max_elem << endl;
  REQUIRE(max_elem == Approx(max_arma));
}


//
