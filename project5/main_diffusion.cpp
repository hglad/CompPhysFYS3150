#include "diffusion.h"

int main(int argc, char const *argv[])
{
  double dx = 0.1;
  double dt = atof(argv[1]);
  double T = 1.;
  int nx = 1./dx;
  int nt = T/dt;

  cout << nt << endl;
  double alpha = dt/(dx*dx);
  cout << alpha << endl;
  double a, c;              // values for sub- and superdiagonal
  //vec d = zeros(nx+2);

  vec b = ones(nx+2)*(1 - 2*alpha);
  vec u = zeros(nx+2);      // current step we want to solve for
  vec y = zeros(nx+2);      // values at previous step

  init_forward(alpha, a, c, b);
//  init_backward(alpha, a, c, b);

  u(0) = 0;   u(nx+1) = 1;
  y(0) = 0;   y(nx+1) = 1;

  // File output
  ofstream myfile;
  string DT = to_string(dt);
  DT = DT.substr(0,6);
  myfile.open ("diffusion_dt=" + DT + ".txt");

  for (int l=1; l < nt; l++)
  {
  //  write_sol(u);
    tridiag(a, c, b, y, u, nx);
    //u(0) = 0;   u(nx+1) = 1;
    /*
    for (int i=1; i < nx-1; i++)
    {
      y(i) = u(i);      // update previous time step
    }
    */
  }
  cout << u << endl;
  return 0;
}
