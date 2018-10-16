#include <iostream>
#include <fstream>
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

// Unit test: test that Earth completes an orbit in a year within a tolerance

// Solvers for planetary orbits in 2D
class solvers
{
  private:
    double pi = M_PI;
    double GM = 4*pi*pi;

  public:
    // Use data matrices to calculate next values in time, given method
    // Verlet
    //
    void Verlet(mat& s, mat& v, mat& a, int i, double M, double dt)
    {

        s(i+1, 0) = s(i, 0) + dt*v(i, 0) + a(i, 0)*dt*dt/2;
        s(i+1, 1) = s(i, 1) + dt*v(i, 1) + a(i, 1)*dt*dt/2;

        double r = sqrt(s(i, 0)*s(i, 0) + s(i, 1)*s(i, 1)); // distance Earth-Sun (should always be 1 AU)
        double F = GM*M/(r*r);    // magnitude of grav. force
        a(i+1, 0) = -F/M * s(i+1, 0)/r;   //cos(theta) = x/r
        a(i+1, 1) = -F/M * s(i+1, 1)/r;

        v(i+1, 0) = v(i, 0) + dt/2*(a(i+1, 0)+a(i, 0));
        v(i+1, 1) = v(i, 1) + dt/2*(a(i+1, 1)+a(i, 1));

        return;
    }
    /*
    void Euler(mat& s, mat& v, mat& a, int i, double M, double dt)
    {
          double r = sqrt(s(i, 0)*s(i, 0) + s(i, 1)*s(i, 1));
          double F = GM*M/(r*r);
          a(i+1, 0) = -F/M * s(i, 0)/r;
          a(i+1, 1) = -F/M * s(i, 1)/r;

          v(i+1, 0) = v(i, 0) + a(i, 0)*dt;
          v(i+1, 1) = v(i, 1) + a(i, 1)*dt;

          s(i+1, 0) = s(i, 0) + v(i, 0)*dt;
          s(i+1, 1) = s(i, 1) + v(i, 1)*dt;

          return;
    }
    */

};
/*
mat Verlet()
{
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
}
*/
