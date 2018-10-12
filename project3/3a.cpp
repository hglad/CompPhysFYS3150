#include <iostream>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  double T = 1;   // simulation time in years
  double dt = 0.001;
  int n = T/dt;

  double pi = atan(1)*4;
  double GM = 4*pi*pi;

  double M_Earth_kg = 6*pow(10, 24);
  double M_Sun_kg = 2*pow(10, 30);
  double M_Sun = 1;
  double M_Earth = M_Earth_kg / M_Sun_kg;

  // Define vectors for physical parameters
  vec x = zeros(n);
  vec y = zeros(n);
  vec v_x = zeros(n);
  vec v_y = zeros(n);
  vec a_x = zeros(n);
  vec a_y = zeros(n);
  x(0) = 1;     // AU
  y(0) = 0;     // AU
  v_y(0) = 2*pi*x(0);   // set initial velocity to 2*pi*r in y-direction

  ofstream myfile;
  char *project3a;
  myfile.open ("project3a.txt");
  for (int i=0; i < n-1; i++)
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
  myfile.close();
  return 0;
}
