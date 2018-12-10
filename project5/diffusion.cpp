#include "diffusion.h"

// Return filename corresponding to method, initialize variables
string init_method(int method, int dim, double dx, double alpha, double &a, double &c, double &b)
{
  string filename;
  string DX = to_string(dx);
  string DIM = to_string(dim);
  DX = DX.substr(0,5);

  if (method == 0)    // Forward
  {
    filename = (DIM + "D_" + "forward_dx=" + DX + ".txt");
  }

  if (method == 1)    // backward
  {
    a = c = -alpha;
    b = (1 + 2*alpha);
    filename = (DIM + "D_" + "backward_dx=" + DX + ".txt");
  }

  if (method == 2)    // crank
  {
    a = c = -alpha;
    b = (2 + 2*alpha);
    filename = (DIM + "D_" + "crank_dx=" + DX + ".txt");
  }
  return filename;
}

void forward(double alpha, vec& u, int nx, double BC1, double BC2)
{
  for (int i=1; i < nx+1; i++)
  {
    u(i) = (1 - 2*alpha)*u(i) + alpha*u(i+1) + alpha*u(i-1);
  }
  set_BCs_1D(u, nx, BC1, BC2);
}

void backward(double a, double c, double b, double alpha, vec& u, vec y, int nx, double BC1, double BC2)
{
  tridiag(a, c, b, y, u, nx);
  set_BCs_1D(u, nx, BC1, BC2);
}

void crank(double a, double c, double b, double alpha, vec& u, vec y, int nx)
{
  for (int i=1; i < nx+1; i++)
  {
    y(i) = alpha*u(i-1) + (2 - 2*alpha)*u(i) + alpha*u(i+1);
  }

  y(0) = 0;
  y(nx+1) = 1;
  //
  tridiag(a, c, b, y, u, nx);

  u(0) = 0;
  u(nx+1) = 1;
}

void forward_2D(double alpha, double D, mat& u, int nx, int ny, double BC1, double BC2)
{
  for (int i=1; i < nx+1; i++)
  {
    for (int j=1; j < ny+1; j++)
    {
      u(i,j) = (1 - 4*D*alpha)*u(i,j) + D*alpha*(u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1));
    }
  }
  set_BCs_2D(u, nx, ny, BC1, BC2);
  return;
}

void forward_source(double fac, double fac2, double dt, mat& u, int nx, int ny, double BC1, double BC2)
{
  double Q;
  for (int i=1; i < nx+1; i++)
  {
    if (i < 20)
    {
      Q = 1.4e-6;
    }

    if ((i <= 40) || (i >= 20))
    {
      Q = 0.35e-6;
    }

    if (i > 40)
    {
      Q = 0.05e-6;
    }

    for (int j=1; j < ny+1; j++)
    {
      u(i,j) = (1 - 4*fac)*u(i,j) + fac*(u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1)) + Q/fac2*dt;
    }
  }
  set_BCs_2D(u, nx, ny, BC1, BC2);
  return;
}

void backward_2D(double a, double c, double b, double alpha, mat& u, mat y, int nx, int ny)
{
  return;
}

void crank_2D(double a, double c, double b, double alpha, mat& u, mat y, int nx, int ny)
{
  return;
}

void set_BCs_1D(vec& u, int nx, double BC1, double BC2)
{
  u(0) = BC1;
  u(nx+1) = BC2;
}

void set_BCs_2D(mat& u, int nx, int ny, double BC1, double BC2)
{
  for (int j=0; j < ny; j++)
  {
    u(0, j) = BC1;          // top
    u(nx+1, j) = BC2;       // bottom
  }
  return;
}

void tridiag(double a, double c, double b, vec y, vec& u, int nx)
{
  // forward substitution
  vec b_tilde = zeros(nx+2); vec y_tilde = zeros(nx+2);
  double a_b_tilde;

  y_tilde(0) = y(0);
  y_tilde(1) = y(1);
  b_tilde(0) = b;
  b_tilde(1) = b;

  for (int i=1; i < nx+1; i++)
  {
    a_b_tilde = a/b_tilde(i-1);	    // save FLOPS by only calculating once
    b_tilde(i) = b - a_b_tilde * c;
    y_tilde(i) = y(i) - a_b_tilde * y_tilde(i-1);
  }

  // backward substitution
  for (int i=nx; i > 0; i--)
  {
    u(i) = (y_tilde(i) - u(i+1)*c)/b_tilde(i);
  }
//  u(nx+1) = (y_tilde(nx) - u(nx+1)*c)/b_tilde(nx+1);
//  y(nx+1) = y(nx+1);
  return;
}

void analytic(int nx, int nt, double L)
{
  double dx = L/(nx+1);
  string DX = to_string(1./nx);
  DX = DX.substr(0,5);
  string filename = ("1D_analytical_dx=" + DX + ".txt");
  ofstream myfile;
  myfile.open(filename);

  double pi = M_PI;
  int N = 1000;

  mat u = zeros(nt, nx+2);
  double An, sum;
  double t_max = 1;
  double t;
  double x;

  double dt = 1./nt;

  for (int l=0; l < nt; l++)
  {
    /*
    u(l, 0) = 0;
    u(l, nx+1) = L;
    */
    t = (double) l/nt;
    for (int i=0; i < nx+2; i++)
    {
      x = i*dx/L;
      sum = 0;

      for (int n=1; n < N+1; n++)
      {
        An = (2/L)*cos(pi*n)/(pi*n*L);

        sum += An * exp(-pow(pi*n/L, 2)*t) * sin(n*pi*x/L);
      }
  //    cout << sum << endl;
      u(l, i) = x + sum;       // i*dx/L: current x-value   u = x + v
  //    cout << u(l, i) << endl;
  //    myfile << l << " " <<  u(i, l) << " ";
    }

    myfile << l << " ";
    for (int i=0; i < nx+2; i++)
    {
      myfile << u(l, i) << " ";
    }
    myfile << endl;


  }

  return;
}
