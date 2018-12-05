#include "diffusion.h"

// Return filename corresponding to method, initialize variables
string init_method(int method, int dim, double dx, double alpha, double &a, double &c, vec &b)
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
    b *= (1 + 2*alpha);
    filename = (DIM + "D_" + "backward_dx=" + DX + ".txt");
  }

  if (method == 2)    // crank
  {
    a = c = -alpha;
    b *= (2 + 2*alpha);
    filename = (DIM + "D_" + "crank_dx=" + DX + ".txt");
  }
  return filename;
}

void forward(double alpha, vec& u, int nx)
{
  for (int i=1; i < nx+1; i++)
  {
    u(i) = (1 - 2*alpha)*u(i) + alpha*u(i+1) + alpha*u(i-1);
  }
  u(0) = 0;
  u(nx+1) = 1;
}

void backward(double a, double c, vec b, double alpha, vec& u, vec y, int nx)
{
  tridiag(a, c, b, y, u, nx);

  u(0) = 0;
  u(nx+1) = 1;
}

void crank(double a, double c, vec b, double alpha, vec& u, vec y, int nx)
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

void forward_2D(double alpha, mat& u, int nx, int ny)
{
  for (int i=1; i < nx+1; i++)
  {
    for (int j=1; j < ny+1; j++)
    {
      u(i,j) = (1 - 4*alpha)*u(i,j) + alpha*(u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1));
    }
  }
  set_BCs_2D(u);
  return;
}

void backward_2D(double a, double c, vec b, double alpha, mat& u, mat y, int nx, int ny)
{
  return;
}

void crank_2D(double a, double c, vec b, double alpha, mat& u, mat y, int nx, int ny)
{
  return;
}

void set_BCs_2D(mat& u)
{
  int nx = u.n_cols;
  int ny = u.n_rows;
  for (int j=0; j < ny; j++)
  {
    u(nx-1, j) = 1;       // bottom
    u(0, j) = 0;          // top
  }
  return;
}

void tridiag(double a, double c, vec b, vec y, vec& u, int nx)
{
  // forward substitution
  vec b_tilde = zeros(nx+2); vec y_tilde = zeros(nx+2);
  double a_b_tilde;

  y_tilde(0) = y(0);
  y_tilde(1) = y(1);
  b_tilde(0) = b(0);
  b_tilde(1) = b(1);

  for (int i=1; i < nx+2; i++)
  {
    a_b_tilde = a/b_tilde(i-1);	    // save FLOPS by only calculating once
    b_tilde(i) = b(i) - a_b_tilde * c;
    y_tilde(i) = y(i) - a_b_tilde * y_tilde(i-1);
  }

  // backward substitution
  for (int i=nx; i > 0; i--)
  {
    u(i) = (y_tilde(i) - u(i+1)*c)/b_tilde(i);
  }
//  u(nx+1) = (y_tilde(nx) - u(nx+1)*c)/b_tilde(nx+1);
//  y(nx+1) = y(nx+1);

}
