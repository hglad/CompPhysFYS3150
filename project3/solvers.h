#ifndef SOLVERS_H
#define SOLVERS_H
#include <armadillo>

using namespace arma;

class Solvers
{
private:
  double pi = M_PI;
  double G = 4*pi*pi;

public:
  vec s;
  vec v;
  vec a;

  double target_M   // mass of body to calculate forces on
  vec M;            // vector of bodies that produce a force on target body

  // Use data matrices to calculate next values in time, given method
  // Verlet
  //
  void Verlet(vec& s, vec& v, vec& a, double M, double dt)
  {

      s(i+1, 0) = s(i, 0) + dt*v(i, 0) + a(i, 0)*dt*dt/2;
      s(i+1, 1) = s(i, 1) + dt*v(i, 1) + a(i, 1)*dt*dt/2;

      double r = sqrt(s(i, 0)*s(i, 0) + s(i, 1)*s(i, 1)); // distance Body-Body (should always be 1 AU)
      double F = G*M/(r*r);    // magnitude of grav. force
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
        double F = G*M/(r*r);
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


#endif /* SOLVERS_H */
