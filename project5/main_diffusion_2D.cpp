#include "diffusion.h"

int main(int argc, char const *argv[])
{
  int method = atoi(argv[1]);

  double h = 0.01;
  double alpha = 0.25;
//  double dt = atof(argv[2]);
  double T = atof(argv[2]);
  double dt = alpha*h*h;
  double L = 1;
//  double alpha = dt/(h*h);
  int nx = L/h;
  int ny = nx;
  int nt = T/dt;
  int add = nt/100;

  cout << nt << endl;
  cout << alpha << endl;
//  double alpha = dt/(dx*dx);
  double a, c;              // values for sub- and superdiagonal
  double b;
  //vec d = zeros(nx+2);

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
  int counter = 0;    // counter used to avoid saving all values
  for (int l=0; l < nt; l++)
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
    if ((l == counter) || (l == nt-1))
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
      counter += 100;   // save every 100th time step
    }


  }
  cout << nt/100 << endl;
//  cout << u << endl;
  myfile.close();
  return 0;
}
