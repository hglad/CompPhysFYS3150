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
  vec d = zeros(nx+2);

  init_forward(alpha, a, c, d);
//  init_backward(alpha, a, c, d);

  //vec d = ones(nx+2)*(1 - 2*alpha);
  vec u = zeros(nx+2);      // current step we want to solve for
  vec y = zeros(nx+2);      // values at previous step

  u(0) = 0;   u(nx+1) = 1;
  y(0) = 0;   y(nx+1) = 1;
  for (int l=1; l < T; l++)
  {
    tridiag(a, c, d, y, u, nx);
    u(0) = 0;   u(nx+1) = 1;
    for (int i=0; i < nx; i++)
    {
      y(i) = u(i);
    }

  }
  cout << u << endl;
  return 0;
}
