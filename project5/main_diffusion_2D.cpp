#include "diffusion.h"

int main(int argc, char const *argv[])
{
  int method = atoi(argv[1]);

  double alpha = 0.25;
  double h = 0.01;
  double dt = alpha*h*h;
  double T = 1; double L = 1;

  int nx = L/h;
  int ny = nx;
  int nt = T/dt;

  cout << nt << endl;
//  double alpha = dt/(dx*dx);
  double a, c;              // values for sub- and superdiagonal
  //vec d = zeros(nx+2);

  vec b = ones(nx+2);
  mat u = zeros(nx+2, ny+2);      // current step we want to solve for
  mat y = zeros(nx+2, ny+2);      // values at previous step

  string filename = init_method(method, 2, h, alpha, a, c, b);

  set_BCs_2D(u);
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
  for (int l=1; l < nt; l++)
  {
    if (method == 0)    // Forward
    {
      forward_2D(alpha, u, nx, ny);
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
    if (l == counter)
    {
      for (int i=0; i < nx+2; i++)
      {
        for (int j=0; j < ny+2; j++)
        {
          myfile << u(i, j) << " ";
        }
        myfile << endl;
      }
      counter += 100;
    }


  }
  cout << u << endl;
  myfile.close();
  return 0;
}
