#include "diffusion.h"

int main(int argc, char const *argv[])
{
  double alpha, D, h, dt, T, L, BC1, BC2;
  double a, c, b;
  int ny, nx, nt, add;
  int method = 0;
  T = 1;
  h = atof(argv[1]);
  alpha = atof(argv[2]);
  //alpha = 0.25;

  BC1 = 0; BC2 = 0;       // top and bottom of grid
//  double dt = atof(argv[2]);
  dt = alpha*h*h;
  cout << dt << endl;
  L = 1;
//  double alpha = dt/(h*h);
  nx = L/h - 2;
  ny = nx;
  nt = T/dt;
  int saved_steps = nt;    // specify number of time steps to write to file
  int save_interval = nt/saved_steps;

  int find_analytic = 1;
  if (find_analytic == 1)
  {
    analytic_2D(nx, ny, nt, L, saved_steps);
  }

  mat u = initial_2D(nx, ny, L);     // current step we want to solve for
  mat y = initial_2D(nx, ny, L);      // values at previous step

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
  int counter = 1;    // counter used to avoid saving all values
  cout << nt << endl;
  for (int l=1  ; l < nt; l++)
  {

    forward_2D(alpha, u, nx, ny, BC1, BC2);
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
