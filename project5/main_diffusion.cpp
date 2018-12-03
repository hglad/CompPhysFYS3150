#include "diffusion.h"

int main(int argc, char const *argv[])
{
  double dx = 0.1;
  double dt = 0.1;

  int nx = 1./dx;
  int nt = 1./dt;
  double T = nt*dt;

  double alpha = dt/(dx*dx);
  double a, c;
  /*
  mat A = init_forward(n, 2);

  double a = A(0,1);
  double d = A(0,0);
  double c = A(1,0);
  */
  a = c = alpha;
  vec d = ones(nx+2)*(1 - 2*alpha);
//  vec d_tilde = d*ones(nx+2);
//  vec f_tilde = ones(nx+2);
  vec u = zeros(nx+2);      // current step we want to solve for
  vec y = zeros(nx+2);      // values at previous step

  u(0) = 0;   u(nx+1) = 1;
  y(0) = 0;   y(nx+1) = 1;
  for (int l=1; l < T; l++)
  {
    tridiag(a, c, d, y, u, nx);
    u(0) = 0;   u(nx+1) = 0;
    for (int i=0; i < nx; i++)
    {
      y(i) = u(i);
    }

  }

  return 0;
}
