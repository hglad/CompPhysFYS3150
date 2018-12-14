#include "diffusion.h"

int main(int argc, char const *argv[])
{
  double dx = atof(argv[1]);
  double alpha = atof(argv[2]);
  double dt = pow(dx,2)*alpha;
  int n = 1./dx - 2;
  int tmax = 1/dt;
  cout << tmax << endl;

  //double alpha = dt/pow(dx,2);
  cout << "Alpha = " << alpha << endl;

  //FESolver(n, alpha, tmax);
  BESolver(n, alpha, tmax, dx, dt);

  //analytical2D(n, dx, 1);
}
