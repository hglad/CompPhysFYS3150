#include "diffusion.h"

int main(int argc, char const *argv[])
{
  double alpha, dx, dt, T, L, BC1, BC2;
  double a, c, b;             // sub-, super- and main diagonal
  int nx, nt, find_analytic, method;
  dx = atof(argv[1]);
  method = atoi(argv[2]);
//  find_analytic = atoi(argv[3]);
  alpha = atof(argv[3]);
  dt = alpha*dx*dx;
//  alpha = 0.5;
//  double dt = 0.001;
  //dt = 0.0001;
  T = 1; L = 1;

  nx = L/dx;
  nt = T/dt;

  int saved_steps = nt;    // specify number of time steps to write to file (subtract 2 because of start and end points)
  int save_interval = (nt/saved_steps);
  cout << save_interval << endl;
//  double alpha = dt/(dx*dx);

  if (find_analytic == 1)
  {
    analytic_1D(nx, nt, L, saved_steps);
  }

  BC1 = 0; BC2 = 1;
  vec u = zeros(nx+2);      // current step we want to solve for
  vec y = zeros(nx+2);      // values at previous step

  string filename = init_method(method, 1, dx, alpha, a, c, b);
  set_BCs_1D(u, nx, BC1, BC2);
  set_BCs_1D(y, nx, BC1, BC2);

  // File output
  ofstream myfile;
  myfile.open(filename);

  for (int i=0; i < nx+2; i++)
  {
    myfile << u(i) << " ";
  }
  myfile << endl;

  // Time loop
  int counter = 1;
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
      crank(a, c, b, alpha, u, y, nx, BC1, BC2);
      y = u;
    }

    // Write results
    if ((l == counter) || (l == nt-1))
    {
      for (int i=0; i < nx+2; i++)
      {
        myfile << u(i) << " ";
      }
      myfile << endl;
      counter += save_interval;
    }
  }
//  cout << u << endl;
  cout << saved_steps << " timesteps saved to file " << filename << endl;
  myfile.close();
  return 0;
}
