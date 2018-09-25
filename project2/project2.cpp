#include "project2.h"

// Function to test diagonalization with armadillo, results can then be compared
// to the results of our Jacobi method
vec arma_diag(mat A)
{
  // Generate vector to hold eigenvalues and matrix to hold eigenvectors
  vec eigval;
  mat eigvec;
  // Decompose into eigenvalues and eigenvalues
  eig_sym(eigval, eigvec, A);
  return sort(eigval);
}

double max_offdiag(mat A, int *k, int *l, double& max_elem)
{
  max_elem = 0;       // reset to 0 so that the if-test below works as intended
  int n = A.n_cols;
  // Find index with highest off-diagonal values, looping through elements that
  // lie to the right of the diagonal
  for (int i=0; i < n; i++)
  {
    for (int j=i+1; j < n; j++)
    {

    if (fabs(A(i, j)) > max_elem)
    {
      max_elem = fabs(A(i, j));
      *k = i;
      *l = j;
    //  printf("%i, %i, %g\n", k, l, max_elem);
    }

    }
  }
  return max_elem;
}

// Function that performs a single Jacobi-rotation
void Jacobi_rot(mat& A, double& max_elem)
{
  double tau, t, t1, t2, s, c;
  int k, l;
  int n = A.n_cols;

  max_elem = max_offdiag(A, &k, &l, max_elem);
  if (A(k,l) != 0.0 )   // Avoid dividing by zero in expression for tau
  {
    tau = (A(l, l) - A(k,k)) / (2*A(k, l));
    t1 = 1./(tau + sqrt(1+tau*tau));
    t2 = -1./(-tau + sqrt(1+tau*tau));
    // Minimize t in order to get smallest possible angle
    if (fabs(t1) <= fabs(t2))
    {
      t = t1;
    }
    else
    {
      t = t2;
    }
    // Calculate sine and cosine from best value of t
    c = 1./sqrt(1 + t*t);
    s = t*c;
  }
  else
  {
    c = 1.0;
    s = 0.0;
  }
  // Generate next matrix by rotation
  double A_kk, A_ll, A_ik, A_il;
  A_kk = A(k, k);   // use A_kk, A_ll as constants to avoid overwriting
  A_ll = A(l, l);

  double cc, ss, cs;
  cc = c*c;    ss = s*s;    cs = c*s;

  A(k,k) = A_kk*cc - 2*A(k, l)*cs + A_ll*ss;
  A(l,l) = A_ll*cc + 2*A(k, l)*cs + A_kk*ss;
  A(k,l) = 0;
  A(l,k) = 0;

  // Loop through elements
  for (int i=0; i < n; i++)
  {
    if (i != k && i != l)
    {
      A_ik = A(i,k);   // define constants to avoid replacing value before it
      A_il = A(i,l);   // is used in the next calculation
      A(i, k) = A_ik*c - A_il*s;
      A(k, i) = A(i, k);        // symmetric matrix
      A(i, l) = A_il*c + A_ik*s;
      A(l, i) = A(i, l);
    //  A(k, l) = (A_kk - A_ll)*c*s + A(k, l)*(c*c - s*s);
    }
  }
  return;
}
