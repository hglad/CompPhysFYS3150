#include "solvers_old.h"
#include "planet_old.h"
#include <fstream>

int main(int argc, char* argv[])
{
  double T = 1;   // simulation time in years
  double dt = 0.001;
  int n = (T/dt) +1; // without +1: planet does not reach starting point
  int np = 1;

  double pi = M_PI;       // define using math library
  double GM = 4*pi*pi;

  double M_Earth_kg = 6*pow(10, 24);
  double M_Sun_kg = 2*pow(10, 30);
  double M_Earth = M_Earth_kg / M_Sun_kg;

  // Create array to hold values for planet at all times
  mat s = zeros(n, 2);
  mat v = zeros(n, 2);
  mat a = zeros(n, 2);

  double M = M_Earth;

  vec s_0 = {1,0};
  // For perfectly circular orbits, velocity should be 2pi*r/T  (AU/yr)
  vec v_0 = {0, 2*pi*s_0(0)/T};

  planet inits;     // class instance for generating initial conditions
  inits.data_mats(n, s, v, a, s_0, v_0, M); // generate initial conds

  ofstream myfile;
  myfile.open ("project3a.txt");
  myfile << s_0(0) << ' ' << s_0(1) << endl;

  solvers Earth_orbit;    // create class instance for solving movement of Earth
  for (int i=0; i < n-1; i++)
  {
    // Compute values at next step i+1 for all planets
    Earth_orbit.Verlet(s, v, a, i, M, dt);
    myfile << s(i+1, 0) << ' ' << s(i+1, 1) << endl;
  }
  myfile.close();

  return 0;
}
