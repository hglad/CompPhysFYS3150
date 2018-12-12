#include "diffusion.h"

int main(int argc, char const *argv[])
{
  double alpha, D, h, dt, T, L, BC1, BC2;
  double a, c, b;
  int ny, nx, nt, add;
  int method = atoi(argv[1]);
  T = atof(argv[2]);
  h = 0.01;
  alpha = 0.25;
  D = 1;                  // diffusion constant
  BC1 = 0; BC2 = 1;       // top and bottom of grid
//  double dt = atof(argv[2]);
  dt = alpha*h*h;
  L = 1;
//  double alpha = dt/(h*h);
  nx = L/h - 2;
  ny = nx;
  nt = T/dt;
  int saved_steps = 500;    // specify number of time steps to write to file
  int save_interval = nt/saved_steps;

  mat u = zeros(nx+2, ny+2);      // current step we want to solve for
  mat y = zeros(nx+2, ny+2);      // values at previous step

  string filename = init_method(method, 2, h, alpha, a, c, b);

  set_BCs_2D(u, nx, ny, BC1, BC2);
//  y(0) = 0;   y(nx+1) = 1;

  // File output
  ofstream myfile;
  myfile.open(filename);

  for (int i=0; i < nx+2; i++)
  {
    for (int j=0; j < ny+2; j++)
    {
      myfile << u(i, j) << " ";
    }
    myfile << endl;
  }

  // Time loop
  int counter = 0;    // counter used to avoid saving all values
  cout << nt << endl;
  for (int l=0; l < nt; l++)
  {
    if (method == 0)    // Forward
    {
      forward_2D(alpha, D, u, nx, ny, BC1, BC2);
    }

    if (method == 1)    // Backward
    {
      backward_2D(a, c, b, alpha, u, y, nx, ny);
      y = u;
    }

    if (method == 2)    // Crank
    {
      crank_2D(a, c, b, alpha, u, y, nx, ny);
      y = u;
    }

    // Write results
    if (l == counter || l == nt-1)
    {
      for (int i=0; i < nx+2; i++)
      {
        for (int j=0; j < ny+2; j++)
        {
          myfile << u(i, j) << " ";
        }
        myfile << endl;
      }
  //    cout << l << "/" << nt << endl;
      counter += save_interval;
    }

  }
  cout << saved_steps << " timesteps saved to file " << filename << endl;
//  cout << u << endl;
  myfile.close();
  return 0;
}
