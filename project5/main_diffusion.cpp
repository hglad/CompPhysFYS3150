#include "diffusion.h"

int main(int argc, char const *argv[])
{
  int method = atoi(argv[1]);

  double alpha = 0.5;
  double dx = 0.1;
  double dt = alpha*dx*dx;
  double T = 1; double L = 1;

  int nx = L/dx;
  int nt = T/dt;

  cout << nt << endl;
//  double alpha = dt/(dx*dx);
  cout << alpha << endl;
  double a, c;              // values for sub- and superdiagonal
  //vec d = zeros(nx+2);

  vec b = ones(nx+2);
  vec u = zeros(nx+2);      // current step we want to solve for
  vec y = zeros(nx+2);      // values at previous step

  string filename = init_method(method, 1, dx, alpha, a, c, b);
//  init_backward(alpha, a, c, b);
//  init_crank(alpha, a, c, b);

  u(0) = 0;   u(nx+1) = 1;
//  y(0) = 0;   y(nx+1) = 1;

  // File output
  analytic(nx, nt);
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
      forward(alpha, u, nx);
    }

    if (method == 1)    // Backward
    {
      backward(a, c, b, alpha, u, y, nx);
    }

    if (method == 2)    // Crank
    {
      crank(a, c, b, alpha, u, y, nx);
    }

    // Write results
    myfile << l << " ";
    for (int i=0; i < nx+2; i++)
    {
      myfile << u(i) << " ";
    }
    myfile << endl;

    y = u;

  }
  cout << u << endl;
  myfile.close();
  return 0;
}
