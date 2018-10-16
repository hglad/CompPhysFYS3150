#include "solvers.h"
#include "planet.h"
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
  double M_Earth = M_Earth_kg / M_Sun_kg; // scale according to solar mass

<<<<<<< HEAD
  // Define vectors for physical parameters
  // np: number of planets
  int np = 1;
  mat s = zeros(n, 2, np);
  mat v = zeros(n, 2, np);
  mat a = zeros(n, 2, np);

  s(0,0,0) = 1;     // AU
  s(0,1,0) = 0;     // AU

  v(0,0,0) = 0;
  // For perfectly circular orbits, velocity should be 2pi*r/T  (AU/yr)
  v(0,1,0) = 2*pi*s(0, 0,0)/T;

  double r_0 = sqrt(s(0,0,0)*s(0,0,0) + s(0,1,0)*s(0,1,0));
  double F_0 = GM*M_Earth/(r_0*r_0);
  a(0,0,0) = -F_0/M_Earth * s(0,0,0)/r_0;
  a(0,1,0) = -F_0/M_Earth * s(0,1,0)/r_0;

  ofstream myfile;
  myfile.open ("project3a.txt");
  myfile << s(0,0,0) << ' ' << s(0,1,0) << endl;

//  planet Earth;           // class instance representing parameters of Earth
=======
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

>>>>>>> 7cab72b862146417a55df1cbe68d5ed74e26ff5b
  solvers Earth_orbit;    // create class instance for solving movement of Earth
  for (int i=0; i < n-1; i++)
  {
<<<<<<< HEAD
    Earth_orbit.Verlet(s, v, a, i, M_Earth, dt);    // Compute values for i+1
    myfile << s(i+1, 0, 0) << ' ' << s(i+1, 1, 0) << endl;
=======
    // Compute values at next step i+1 for all planets
    Earth_orbit.Verlet(s, v, a, i, M, dt);
    myfile << s(i+1, 0) << ' ' << s(i+1, 1) << endl;
>>>>>>> 7cab72b862146417a55df1cbe68d5ed74e26ff5b
  }
  myfile.close();

  return 0;




}
