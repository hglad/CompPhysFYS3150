#include "diffusion.h"

// Return filename corresponding to method, initialize variables
string init_method(int method, int dim, double dx, double alpha, double &a, double &c, double &b)
{
  string filename;
  string DX = to_string(dx);
  string DIM = to_string(dim);
  string ALPHA = to_string(alpha);
  ALPHA = ALPHA.substr(0,4);
  DX = DX.substr(0,5);

  if (method == 0)    // Forward
  {
    filename = (DIM + "D_" + "forward_dx=" + DX + "_alpha=" + ALPHA + ".txt");
  }

  if (method == 1)    // backward
  {
    a = c = -alpha;
    b = (1 + 2*alpha);
    filename = (DIM + "D_" + "backward_dx=" + DX + "_alpha=" + ALPHA + ".txt");
  }

  if (method == 2)    // crank
  {
    a = c = -alpha;
    b = (2 + 2*alpha);
    filename = (DIM + "D_" + "crank_dx=" + DX + "_alpha=" + ALPHA + ".txt");
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
  set_BCs_1D(y, nx, BC1, BC2);
  tridiag(a, c, b, y, u, nx);
  set_BCs_1D(u, nx, BC1, BC2);
}

void crank(double a, double c, double b, double alpha, vec& u, vec y, int nx, double BC1, double BC2)
{
  for (int i=1; i < nx+1; i++)
  {
    y(i) = alpha*u(i-1) + (2 - 2*alpha)*u(i) + alpha*u(i+1);
  }

  set_BCs_1D(y, nx, BC1, BC2);
  //
  tridiag(a, c, b, y, u, nx);

  set_BCs_1D(u, nx, BC1, BC2);
}

void forward_2D(double alpha, mat& u, int nx, int ny, double BC1, double BC2)
{
  for (int i=1; i < nx+1; i++)
  {
    for (int j=1; j < ny+1; j++)
    {
      u(i,j) = (1 - 4*alpha)*u(i,j) + alpha*(u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1));
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
  for (int j=0; j < ny+2; j++)
  {
    u(0, j) = BC1;          // top
    u(nx+1, j) = BC2;       // bottom
    u(j, 0) = BC1;          // top
    u(j, ny+1) = BC2;       // bottom
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
  return;
}

mat initial_2D(int nx, int ny, double L)
{
  mat u = zeros(nx+2, ny+2);
  double x, y, dx, dy;
  double pi = M_PI;
  dx = L/(nx+1);
  dy = dx;
  for (int i=0; i < nx+2; i++)
  {
    x = i*dx;
    for (int j=0; j < ny+2; j++)
    {
      y = j*dy;
      u(i,j) = sin(pi*x) * sin(pi*y);
    }
  }
  return u;
}

void analytic_1D(int nx, int nt, double L, int saved_steps)
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
  u(0, nx+1) = 1;
  double An, sum;
  double t_max = 1;
  double t;
  double x;

  double dt = 1./nt;
  int save_interval = nt/saved_steps;

  myfile << 0 << " ";
  for (int i=0; i < nx+2; i++)
  {
    myfile << u(0, i) << " ";
  }
  myfile << endl;

  int counter = 1;
  for (int l=1; l < nt; l++)
  {

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
      u(l, i) = x + sum;       // i*dx/L: current x-value   u = x + v

    }
    if ((l == counter) || (l == nt-1))
    {
      myfile << l << " ";
      for (int i=0; i < nx+2; i++)
      {
        myfile << u(l, i) << " ";
      }
      myfile << endl;
      counter += save_interval;

    }

  }

  return;
}

void analytic_2D(int nx, int ny, int nt, double L, int saved_steps)
{
  double dx = L/(nx+1);
//  double dx = h;
  double dy = dx;
  double BC1, BC2;
  BC1 = BC2 = 0;

  string DX = to_string(dx);
  DX = DX.substr(0,5);
  string filename = ("2D_analytical_dx=" + DX + ".txt");
  ofstream myfile;
  myfile.open(filename);

  double pi = M_PI;

  mat u = initial_2D(nx, ny, L);
  //BC1 = 0; BC2 = 1;
  set_BCs_2D(u, nx, ny, BC1, BC2);

  double t_max = 1;
  double t;
  double x, y;

  double dt = 1./nt;
  int save_interval = nt/saved_steps;

  for (int i=0; i < nx+2; i++)
  {
    for (int j=0; j < ny+2; j++)
    {
      myfile << u(i, j) << " ";
    }
    myfile << endl;
  }

  int counter = 1;
  int write;
  for (int l=1; l < nt; l++)
  {
    write = 0;
    if (l == counter || l == nt-1)
    {
      write = 1;
    }

    t = (double) l/nt;

    for (int i=1; i < nx+1; i++)
    {
      x = i*dx;

      for (int j=1; j < ny+1; j++)
      {
        y = j*dy;

        u(i,j) = sin(pi*x) * sin(pi*y) * exp(-2*pi*pi*t);
      }
    }
    //set_BCs_2D(u, nx, ny, BC1, BC2);
    if (write == 1)
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

  return;
}

void JSolver(double e, double d, int n, double alpha, mat &u, mat rhs){
        int maxIter = 1000;
        double diff = 1.0;
        int iter = 0;
        double tol = 1E-10;
        mat u_old;

        while( iter < maxIter && diff > tol) {
                diff = 0;
                u_old = u;
                for(int i = 1; i < n+1; i++) {
                        for(int j = 1; j < n+1; j++) {
                                u(i,j) = 1.0/d*(alpha*(u_old(i+1,j) + u_old(i,j+1) +
                                                       u_old(i-1,j)+ u_old(i,j-1)) + rhs(i,j));
                                diff += abs(u(i,j) - u_old(i,j));
                        }
                }
                iter++;
        }
  //      cout << iter << endl;
}

void BESolver(int n, double alpha, int tmax, double dx, double dt){
        mat u = initial_2D(n, n, 1);  // Au = r
        double BC1, BC2;
        BC1 = BC2 = 0;
        // Boundary conditions (u(0) set by zeros)
        set_BCs_2D(u, n, n, BC1, BC2);

        // Matrix elements of tridiagonal matrix
        double e = -alpha;
        double d = 1.0 + 4.0*alpha;

        int saved_steps = tmax;    // specify number of time steps to write to file
        int save_interval = tmax/saved_steps;

        mat A = zeros<mat>(n+2,n+2);
        A.diag() += d; // center diagonal
        A.diag(1) += e; // upper diagonal
        A.diag(-1) += e; // lower diagonal

        string DX = to_string(dx);
        DX = DX.substr(0,5);
        string ALPHA = to_string(alpha);
        ALPHA = ALPHA.substr(0,4);

        string filename = ("2D_backward_dx=" + DX + "_alpha=" + ALPHA + ".txt");

        ofstream outfile;
        outfile.open(filename);

        // writing initial state to file
        for(int i = 0; i < n+2; i++) {
                for(int j = 0; j < n+2; j++) {
                        outfile << u(i,j) << " ";
                }
                outfile << endl;
        }

        int counter = 1;
        cout << tmax << endl;
        for(int j = 1; j < tmax; j++) {
                JSolver(e, d, n, alpha, u, u);

                // Preserving boundary conditions
                set_BCs_2D(u, n, n, BC1, BC2);

                // writing to file
                if(j == counter || j == tmax-1) {
                        for(int i = 0; i < n+2; i++) {
                                for(int j = 0; j < n+2; j++) {
                                        outfile << u(i,j) << " ";
                                }
                                outfile << endl;
                        }
                        counter += save_interval;
                }
        }
        outfile.close();

}
