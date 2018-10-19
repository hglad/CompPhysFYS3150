#ifndef SOLVERS_H
#define SOLVERS_H

#include "project3.h"
#include "Body.h"

class Solvers
{
private:
  double pi = M_PI;
  double G = 4*pi*pi;

public:
  Body target_Body;   // object of body to calculate forces on
  vector<Body> Bodies;          // vector of body objects that produce a force on target body
  vec old_Pos;
  vec old_Vel;
  vec acc = {0,0};  // initialize acceleration as zero, add contributions later

  double delta_x;
  double delta_y;
  double old_d;
  double temp_d;
  double temp_a;
  int n_ext_bodies;

  double dt;

  Solvers();

  Solvers(Body target_body, vector<Body> bodies, double DT);

  void Verlet(vec& pos, vec& vel, vec& old_Acc);

};

#endif /* SOLVERS_H */

  /*
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
