// Compile: g++ project2.cpp test-project2b.cpp -o test.x -DARMA_DONT_USE_WRAPPER -lblas -llapack
// Test that our method finds the eigenvalues of a random 5x5 matrix
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "project2.h"

vec jacobi_diag(mat& B)
{
  int m = B.n_cols;     // dimension of matrix
  double eps = 1e-8;
  double max_elem = 2*eps; // initialize max variable with arbitrary value > eps

  while (max_elem > eps)
  {
    Jacobi_rot(B, max_elem);
  }
  // Store eigenvalues in vector
  vec eigvals = B.diag();

  cout << sort(eigvals) << endl;

  // Return eigenvalues sorted from lowest to highest in value (like Armadillo does)
  return sort(eigvals);
}

// Test for a random tridiagonal matrix. Dense matrix produces less accurate results.
mat initialize(int m)
{
  srand (time(NULL));
  double d = rand() % 1000;
  double a = rand() % 1000;
  mat B = generate_mat(m, d, a, -a);
  return B;
}

int m = 5;
mat B = initialize(m);

// Test our results
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
  B.diag().fill(0);   // do not check main diagonal for elements
  B = abs(B);
  int k, l;
  double max_arma = B.max();
  double max_elem = max_offdiag(B, &k, &l, max_elem);

  cout << max_arma << ' ' << max_elem << endl;
  REQUIRE( abs(max_elem) == Approx(max_arma) );
}

//
