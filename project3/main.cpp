#include "Body.h"
#include "math.h"
#include "Solvers.h"
#include <vector>

int main(int argc, char const *argv[])
{
  double pi = M_PI;
  double M_Sun_kg = 1.989*pow(10,30);
  double MSun = 1;

  double ME = 5.972*pow(10,24)/M_Sun_kg;

  vec PosSun = {0,0};
  vec VelSun = {0,0};

  vec PosE = {1,0};
  vec VelE = {0, 2*pi};

  Body Sun(PosSun, VelSun, MSun);
  Body Earth(PosE, VelE, ME);

  cout << Earth.distance(Sun) << endl;

  double T = 1;      //years
  double dt = 0.01;   //years
  int n = (T/dt)+1;

  double dOld;
  double d;
  double FOld;
  double F;

  vec PosOld;
  vec VelOld;

  vec AccE = {0,0};
  vec AccOld;
  cout << n << endl;

  ofstream myfile;
  myfile.open ("project3.txt");
  myfile << PosE(0) << ' ' << PosE(1) << endl;

  vector<Body> bodies = {Sun};    // vector with elements that are Body objects

  vec old_Acc = {0,0};
  Solvers Earth_Sun(Earth, bodies, dt);

  // Integration loop
  for (int i=0; i < n; i++)
  {
    Earth_Sun.Verlet(PosE, VelE, old_Acc);
  }

  //cout << Earth.distance(Sun) << endl;
  //cout << PosE(0) << ' ' << PosE(1) << endl;
  myfile.close();

  // Verlet
  /*
  for (int i=0; i < n; i++)
  {
    PosOld = PosE;
    VelOld = VelE;

    dOld = Earth.distance(Sun);       // compute distance Earth-Sun
    FOld = Earth.grav_force(Sun);
    AccOld = {-FOld/ME * PosOld(0)/dOld, -FOld/ME * PosOld(1)/dOld};

    PosE(0) = PosOld(0) + dt*VelOld(0) + AccOld(0)*dt*dt/2;
    PosE(1) = PosOld(1) + dt*VelOld(1) + AccOld(1)*dt*dt/2;

    F = Earth.grav_force(Sun); d = Earth.distance(Sun);

    AccE(0) = -F/ME * PosE(0)/d;
    AccE(1) = -F/ME * PosE(1)/d;

    VelE(0) = VelOld(0) + dt/2*( AccE(0) + AccOld(0) );
    VelE(1) = VelOld(1) + dt/2*( AccE(1) + AccOld(1) );
    myfile << PosE(0) << ' ' << PosE(1) << endl;

  }
  */
  return 0;

}
