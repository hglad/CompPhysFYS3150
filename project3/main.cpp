#include "Body.h"
#include "math.h"
#include "Solvers.h"
#include <vector>

int main(int argc, char const *argv[])
{
  // Simulation parameters
  double T = 1;      //years
  double dt = 0.01;   //years
  int n = (T/dt)+1;
  cout << n << endl;

  // Physical constants
  double pi = M_PI;
  double M_Sun_kg = 1.989*pow(10,30);
  double MSun = 1;

  double ME = 5.972*pow(10,24)/M_Sun_kg;

  // Initial vectors
  vec PosSun = {0,0};
  vec VelSun = {0,0};

  vec PosE = {1,0};
  vec VelE = {0, 2*pi};

  // Create objects representing planets/sun
  Body Sun(PosSun, VelSun, MSun);
  Body Earth(PosE, VelE, ME);

  cout << Earth.distance(Sun) << endl;

  ofstream myfile;
  myfile.open ("project3.txt");
  myfile << PosE(0) << ' ' << PosE(1) << endl;

  vector<Body> bodies = {Sun};    // vector with elements that are Body objects

  Solvers Earth_Sun(Earth, bodies, dt);

  // Integration loop
  vec old_Acc = {0,0};
//  double a_0 = Earth.grav_force(Sun)/ME;
//  old_Acc = {a_0*1, a_0*0};
  for (int i=0; i < n; i++)
  {
    Earth_Sun.Verlet(PosE, VelE, old_Acc);
    //Earth_Sun.Euler(PosE, VelE);
    myfile << PosE(0) << ' ' << PosE(1) << endl;
  //  myfile << Earth.Pos(0) << ' ' << Earth.Pos(1) << endl;
    // Update object values
    Earth.update(PosE, VelE);
  }

  cout << Earth.Pos << endl;
  myfile.close();

  return 0;

}
