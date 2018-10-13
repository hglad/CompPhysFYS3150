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

  int n = (T/dt) +1; // without +1: planet does not reach starting point
  cout << n << endl;

  double pi = M_PI;       // define using math library
  double GM = 4*pi*pi;

  double M_Earth_kg = 6*pow(10, 24);
  double M_Sun_kg = 2*pow(10, 30);
  double M_Sun = 1;       // use solar mass as scaling factor
  double M_Earth = M_Earth_kg / M_Sun_kg;

  // Define vectors for physical parameters
  mat s = zeros(n, 2);
  mat v = zeros(n, 2);
  mat a = zeros(n, 2);

  s(0,0) = 1;     // AU
  s(0,1) = 0;     // AU

  v(0,0) = 0;
  v(0,1) = 2*pi*s(0, 0);   // set initial velocity to 2*pi*r in y-direction

  double r_0 = sqrt(s(0,0)*s(0,0) + s(0,1)*s(0,1));
  double F_0 = GM*M_Earth/(r_0*r_0);
  a(0,0) = -F_0/M_Earth * s(0,0);
  a(0,1) = -F_0/M_Earth * s(0,1);

  ofstream myfile;
  myfile.open ("project3a.txt");
  myfile << s(0,0) << ' ' << s(0,1) << endl;
  // Euler
  /*
  for (int i=0; i < n-1; i++)
  {
      double r = sqrt(x(i)*x(i) + y(i)*yaa(i));
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
  for (int i=0; i < n-1; i++)
  {
      s(i+1, 0) = s(i, 0) + dt*v(i, 0) + a(i, 0)*dt*dt/2;
      s(i+1, 1) = s(i, 1) + dt*v(i, 1) + a(i, 1)*dt*dt/2;

      double r = sqrt(s(i, 0)*s(i, 0) + s(i, 1)*s(i, 1)); // distance Earth-Sun (should always be 1 AU)
      double F = GM*M_Earth/(r*r);    // magnitude of grav. force
      a(i+1, 0) = -F/M_Earth * s(i+1, 0);     // cos(theta) = x if radius is 1 AU
      a(i+1, 1) = -F/M_Earth * s(i+1, 1);

      v(i+1, 0) = v(i, 0) + dt/2*(a(i+1, 0)+a(i, 0));
      v(i+1, 1) = v(i, 1) + dt/2*(a(i+1, 1)+a(i, 1));

      myfile << s(i+1, 0) << ' ' << s(i+1, 1) << endl;
  }

  myfile.close();

  return 0;
}
