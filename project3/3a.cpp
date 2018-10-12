#include <iostream>
#include <fstream>
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  double T = 1;   // simulation time in years
  double dt = 0.001;
  int n = (T/dt);
  cout << n << endl;

  double pi = M_PI;       // from math.h
  double GM = 4*pi*pi;

  double M_Earth_kg = 6*pow(10, 24);
  double M_Sun_kg = 2*pow(10, 30);
  double M_Sun = 1;
  double M_Earth = M_Earth_kg / M_Sun_kg;

  // Define vectors for physical parameters
  // with zeros(n): planet does not reach starting point
  vec x = zeros(n+1);
  vec y = zeros(n+1);
  vec v_x = zeros(n+1);
  vec v_y = zeros(n+1);
  vec a_x = zeros(n+1);
  vec a_y = zeros(n+1);
  x(0) = 1;     // AU
  y(0) = 0;     // AU
  v_y(0) = 2*pi*x(0);   // set initial velocity to 2*pi*r in y-direction

  double F_0 = GM*M_Earth/(x(0)*x(0));
  a_x(0) = -F_0/M_Earth * x(0);
  a_y(0) = -F_0/M_Earth * y(0);

  ofstream myfile;
  char *project3a;
  myfile.open ("project3a.txt");
  myfile << x(0) << ' ' << y(0) << endl;
  // Euler
  /*
  for (int i=0; i < n; i++)
  {
      double r = sqrt(x(i)*x(i) + y(i)*y(i));
      double F = GM*M_Earth/(r*r);
      a_x(i+1) = -F/M_Earth * x(i);
      a_y(i+1) = -F/M_Earth * y(i);

      v_x(i+1) = v_x(i) + a_x(i+1)*dt;
      v_y(i+1) = v_y(i) + a_y(i+1)*dt;

      x(i+1) = x(i) + v_x(i+1)*dt;
      y(i+1) = y(i) + v_y(i+1)*dt;
      myfile << x(i) << ' ' << y(i) << endl;
  }
  */
  // Verlet
  for (int i=0; i < n; i++)
  {
      x(i+1) = x(i) + dt*v_x(i) + a_x(i)*dt*dt/2;
      y(i+1) = y(i) + dt*v_y(i) + a_y(i)*dt*dt/2;

      double r = sqrt(x(i)*x(i) + y(i)*y(i)); // distance Earth-Sun (should always be 1 AU)
      double F = GM*M_Earth/(r*r);    // magnitude of grav. force
      a_x(i+1) = -F/M_Earth * x(i+1);
      a_y(i+1) = -F/M_Earth * y(i+1);

      v_x(i+1) = v_x(i) + dt/2*(a_x(i+1)+a_x(i));
      v_y(i+1) = v_y(i) + dt/2*(a_y(i+1)+a_y(i));

      myfile << x(i+1) << ' ' << y(i+1) << endl;
  }

  myfile.close();

  return 0;
}
