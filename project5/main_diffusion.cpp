#include "diffusion.h"

int main(int argc, char const *argv[])
{
  double alpha, dx, dt, T, L, BC1, BC2;
  double a, c, b;             // sub-, super- and main diagonal
  int nx, nt;
  int method = atoi(argv[1]);

  alpha = 0.5;
//  double dx = 1./120.;
  dx = 0.1;
//  double dt = 0.001;
  dt = alpha*dx*dx;
  //dt = 0.0001;
  T = 1; L = 1;

  nx = L/dx;
  nt = T/dt;
//    double alpha = dt/(dx*dx);
  analytic(nx, nt, L);

  BC1 = 0; BC2 = 1;
  vec u = zeros(nx+2);      // current step we want to solve for
  vec y = zeros(nx+2);      // values at previous step

  string filename = init_method(method, 1, dx, alpha, a, c, b);
  set_BCs_1D(u, nx, BC1, BC2);
  set_BCs_1D(y, nx, BC1, BC2);

  // File output
  ofstream myfile;
  myfile.open(filename);

  myfile << 0 << " ";
  for (int i=0; i < nx+2; i++)
  {
    myfile << u(i) << " ";
  }
  myfile << endl;

  // Time loop
  for (int l=1; l < nt; l++)
  {
    if (method == 0)    // Forward
    {
      forward(alpha, u, nx, BC1, BC2);
    }

    if (method == 1)    // Backward
    {
      backward(a, c, b, alpha, u, y, nx, BC1, BC2);
      y = u;
    }

    if (method == 2)    // Crank
    {
      crank(a, c, b, alpha, u, y, nx);
      y = u;
    }

    // Write results
    myfile << l << " ";
    for (int i=0; i < nx+2; i++)
    {
      myfile << u(i) << " ";
    }
    myfile << endl;

  }
  cout << u << endl;
  myfile.close();
  return 0;
}
