#include "solvers.h"
#include "planet.h"

int main(int argc, char* argv[])
{
  double T = 1;   // simulation time in years
  double dt = 0.001;

  int n = (T/dt) +1; // without +1: planet does not reach starting point

  double pi = M_PI;       // define using math library
  double GM = 4*pi*pi;

  double M_Earth_kg = 6*pow(10, 24);
  double M_Sun_kg = 2*pow(10, 30);
  double M_Earth = M_Earth_kg / M_Sun_kg;

  // Define vectors for physical parameters
  mat s = zeros(n, 2);
  mat v = zeros(n, 2);
  mat a = zeros(n, 2);

  s(0,0) = 1;     // AU
  s(0,1) = 0;     // AU

  v(0,0) = 0;
  // For perfectly circular orbits, velocity should be 2pi*r/T  (AU/yr)
  v(0,1) = 2*pi*s(0, 0)/T;

  double r_0 = sqrt(s(0,0)*s(0,0) + s(0,1)*s(0,1));
  double F_0 = GM*M_Earth/(r_0*r_0);
  a(0,0) = -F_0/M_Earth * s(0,0)/r_0;
  a(0,1) = -F_0/M_Earth * s(0,1)/r_0;

  ofstream myfile;
  myfile.open ("project3a.txt");
  myfile << s(0,0) << ' ' << s(0,1) << endl;

  planet Earth;           // class instance representing parameters of Earth (not used yet)
  solvers Earth_orbit;    // create class instance for solving movement of Earth

  for (int i=0; i < n-1; i++)
  {
    Earth_orbit.Verlet(s, v, a, i, M_Earth, dt);    // Compute values at i+1
    myfile << s(i+1, 0) << ' ' << s(i+1, 1) << endl;
  }
  myfile.close();

  return 0;
}
