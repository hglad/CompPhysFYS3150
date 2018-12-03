#include "diffusion.h"

mat init_forward(int n, double alpha)
{
  mat A = zeros(n,n);

  for (int i=0; i < n; i++)
  {
    A(i, i)    = 1 - 2*alpha;      // Main diagonal
    A(i, i+1)  = alpha;      // Upper diagonal values
    A(i+1, i)  = alpha;      // Lower diagonal values
  }
  return A;
}

void tridiag(double a, double c, vec d, vec& f, vec& u, int n)
{
  // forward
  for (int i=1; i < n-1; i++)
  {
    d(i) = d(i) - a*c/d(i-1);
    f(i) = f(i) - a/d(i-1)*f(i-1);
  }
  // backward
  for (int i=n-1; i > 1; i--)
  {
    u(i) = (f(i) - u(i+1)*c)/d(i);
  }
}
