#include "diffusion.h"

int main(int argc, char const *argv[])
{
  int method = 0;
  double L = 120*1000;
  double h = L/120;
  double alpha = 0.25;
//  double dt = atof(argv[2]);
  double T = atof(argv[1]);
  double dt = alpha*h*h;

//  double alpha = dt/(h*h);
  int nx = L/h;
  int ny = nx;
  int nt = 1e5;
  int add = nt/100;

  cout << nt << endl;
  cout << alpha << endl;
//  double alpha = dt/(dx*dx);
  double a, c;              // values for sub- and superdiagonal
  //vec d = zeros(nx+2);

  vec b = ones(nx+2);
  mat u = zeros(nx+2, ny+2);      // current step we want to solve for
  mat y = zeros(nx+2, ny+2);      // values at previous step

  string DX = to_string(h);
  DX = DX.substr(0,5);

  string filename = ("lithos_forward_dx=" + DX + ".txt");
  set_BCs_source(u);
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

  // Thermal constants - non-scaled
  double k = 2.5;
  double rho = pow(3.510, 3);
  double Cp = 1000;
  double D = k/(rho*Cp);
  double fac = alpha*D;
  double fac2 = rho*Cp;

  // Time loop
  int counter = 0;    // counter used to avoid saving all values
  for (int l=0; l < nt; l++)
  {
    forward_source(fac, fac2, dt, u, nx, ny);

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
      cout << l << "/" << nt << endl;
      counter += 5000;   // save every 100th time step
    }

  }
//  cout << u << endl;
  myfile.close();
  return 0;
}
