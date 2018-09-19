#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

// Function to test diagonalization with armadillo
void arma_diag_test(mat A)
{
  // Generate vector to hold eigenvalues and matrix to hold eigenvectors
  vec eigval;
  mat eigvec;
  // Decompose into eigenvalues and eigenvalues
  eig_sym(eigval, eigvec, A);
  cout << "Eigenvalues and eigenvectors:" << endl;
  cout << eigval << eigvec << endl;
  return;
}

// Calculate sine and cosine values, given matrix and indices to max value
double trig(mat& A, double k, double l)
{
  double tau, t, t1, t2, s, c;
  tau = (A(l, l) - A(k,k)) / 2*A(k, l);
  t1 = -tau + sqrt(1 + tau*tau);
  t2 = -tau - sqrt(1 + tau*tau);

  // Minimize t in order to get smallest possible angle
  if (fabs(t1) <= fabs(t2))
  {
    t = t1;
  }
  else
  {
    t = t2;
  }

  c = 1./sqrt(1 + t*t);
  s = t*c;
  return s, c;
}
// Perform one Jacobi rotation
void Jacobi_rot(mat& A, int n, double& max_elem)
{
  max_elem = 0; // if max_elem is large enough, set to 0 before finding next max
  double A_kk, A_ll, A_ik, A_il;
  int k, l;   // indices for maximum off-diagonal values

  // Find index with highest off-diagonal values
  for (int i=0; i < n-1; i++)
  {
    for (int j=0; i+1 < n-1; j++)
    {

    if (fabs(A(i, j)) > max_elem)
    {
      max_elem = fabs(A(i, j));
      k = i;
      l = j;
    //  printf("max_elem in rot-function: %g, i=%i, j=%i\n", max_elem, k, l);
    }

    }
  }

  // Calculate sine and cosine with own function
  double s, c = trig(A, k, l);

  // Generate next matrix by rotation
  A_kk = A(k, k);   // use A_kk, A_ll as constants to avoid overwriting
  A_ll = A(l, l);
  A(k,l) = 0;
  A(l,k) = 0;

  A(k,k) = A_kk*c*c - 2*A(k, l)*c*s + A_ll*s*s;
  A(l,l) = A_ll*c*c + 2*A(k, l)*c*s + A_kk*s*s;
//  A(k,l) = (A_kk - A_ll)*c*s + A(k, l)*(c*c - s*s);

  // Loop through elements
  for (int i=0; i < n; i++)
  {
    if (i != k && i != l)
    {
      A_ik = A(i,k);   // define constants to avoid replacing value in algo
      A_il = A(i,l);
      A(i, k) = A_ik*c - A_il*s;
      A(i, l) = A_il*c + A_ik*s;
      A(l, l) = A_ll*c*c + 2*A(k,l)*c*s + A_kk*s*s;
    //  A(k, l) = (A_kk - A_ll)*c*s + A(k, l)*(c*c - s*s);
    //  A(k, i) = A(i, k);
    //  A(l, i) = A(i, l);
    }
  }
  return;
}

int main(int argc, char *argv[])
{
  // Initialize matrix and necessary variables
  int n = atoi(argv[1]);
  double eps = 1e-8;    // tolerance to represent 0

  mat A = zeros(n,n);
  double h = 1./(n+1);
  double d = 2./(h*h);
  double a = -1./(h*h);

  // Generate tridiagonal matrix
  for (int i=0; i < n-1; i++)
  {
    A(i, i)    = d;      // Main diagonal
    A(i, i+1)  = a;      // Upper diagonal values
    A(i+1, i)  = a;      // Lower diagonal values
  }
  A(n-1, n-1) = d;

  arma_diag_test(A);

  // Rotate matrix until value of highest element is close to 0
  cout << A << endl;    // print initial matrix
  double max_elem = 2*eps;  // initialize max variable with arbitrary value > eps
  int num_rotations = 0;
  while (max_elem > eps)
  {
    Jacobi_rot(A, n, max_elem);
    num_rotations += 1;
  //  printf("max_elem outside function: %g\n", max_elem);
  }

  cout << A << endl;    // print resulting matrix
  printf("Number of rotations performed: %i\n", num_rotations );
  return 0;

}
