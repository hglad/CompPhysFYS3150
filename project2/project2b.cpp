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
  eig_sym(eigval, eigvec, A);
  cout << "Eigenvalues and eigenvectors:" << endl;
  cout << eigval << eigvec << endl;
  return;
}

void Jacobi_rot(mat& A, int n)
{

  double tau, t, c, s, A_kk, A_ll, A_ik, A_il;
  //Jacobi's method
  double eps = 1e-8;
  double max_elem = 0;
  int k, l;   // indices for maximum off-diagonal values

  // Find index with highest off-diagonal values
  for (int i=0; i < n-1; i++)
  {
    for (int j=0; j < n-1; j++)
    {

    if (fabs(A(i, j+1)) > max_elem)
    {
      max_elem = A(i, j+1);
      k = i;
      l = j+1;
    }
    //cout << max_elem << endl;

    }
  }
  // Calculate values needed for rotation
  double t1, t2;
  tau = (A(l, l) - A(k,k)) / 2*A(k, l);
  t1 = -tau + sqrt(1 + tau*tau);
  t2 = -tau - sqrt(1 + tau*tau);

  // Minimize t in order to get smallest possible angle
  if (abs(t1) <= abs(t2))
  {
    t = t1;
  }
  else
  {
    t = t2;
  }

  c = 1./sqrt(1 + t*t);
  s = t*c;

  // Generate next matrix by rotation
  A_kk = A(k, k);   // use A_kk, A_ll as constants to avoid overwriting
  A_ll = A(l, l);

  A(k,k) = A_kk*c*c - 2*A(k, l)*c*s + A(l, l)*s*s;
  A(l,l) = A_ll*c*c + 2*A(k, l)*c*s + A(k, k)*s*s;

//  A(k,l) = (A_kk - A_ll)*c*s + A(k, l)*(c*c - s*s);
  A(k,l) = 0;
  A(l,k) = 0;

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
//  int m = atoi(argv[2]);    // number of times to perform rotation
  int m = 5*n*n;
  cout << m << endl;

  mat A = zeros(n,n);
  double h = 1./(n+1);
  double d = 2./(h*h);
  double a = -1./(h*h);

  for (int i=0; i < n-1; i++)
  {
    A(i, i)    = d;      // Main diagonal
    A(i, i+1)  = a;      // Upper diagonal values
    A(i+1, i)  = a;      // Lower diagonal values
  }
  A(n-1, n-1) = d;

  arma_diag_test(A);

  // Rotate matrix m times
  cout << A << endl;    // print initial matrix
  for (int j=0; j < m; j++)
  {
    Jacobi_rot(A, n);
  }

  cout << A << endl;    // print resulting matrix

  return 0;

}
