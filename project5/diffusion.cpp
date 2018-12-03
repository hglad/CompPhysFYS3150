#include "diffusion.h"

void init_forward(double alpha, double &a, double &c, vec& d)
{
  a = c = alpha;
  d *= (1 - 2*alpha);

  return;
}

void init_backward(double alpha, double &a, double &c, vec& d)
{
  a = c = -alpha;
  d *= (1 + 2*alpha);

  return;
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
