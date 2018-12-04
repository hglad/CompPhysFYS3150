#include "diffusion.h"

void init_forward(double alpha, double &a, double &c, vec& b)
{
  a = c = alpha;
  //b *= (1 - 2*alpha);
  return;
}

void init_backward(double alpha, double &a, double &c, vec& b)
{
  a = c = -alpha;
  //b *= (1 - 2*alpha);
  return;
}

void tridiag(double a, double c, vec b, vec& y, vec& u, int nx)
{
  // forward substitution
  vec b_tilde = zeros(nx+2); vec y_tilde = zeros(nx+2);
  double a_b_tilde;

  y_tilde(0) = y(0);
  y_tilde(1) = y(1);
  b_tilde(0) = b(0);
  b_tilde(1) = b(1);

  for (int i=2; i < nx+2; i++)
  {
    a_b_tilde = a/b_tilde(i-1);	    // save FLOPS by only calculating once
    b_tilde(i) = b(i) - a_b_tilde * c;
    y_tilde(i) = y(i) - a_b_tilde * y_tilde(i-1);
  //  b(i) = b(i) - a*c/b(i-1);           // main diagonal
  //  y(i) = y(i) - a/b(i-1)*y(i-1);      // RHS
  }

  // backward substitution
  for (int i=nx; i > 0; i--)
  {
  //  u(i) = (y(i) - u(i+1)*c)/b(i);
    u(i) = (y_tilde(i) - u(i+1)*c)/b_tilde(i);
  //  y(i) = u(i);      // update RHS
  }

}
/*
void write_sol(vec u, string filename)
{

}
*/
