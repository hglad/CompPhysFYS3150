#include "Body.h"
#include <armadillo>

// define constructors for Body
Body::Body()   // generate an Earth body 1 AU from the Sun with circular orbit
{
  double M_Sun_kg = 1.989*pow(10,30);
  double ME = 5.972*pow(10,24)/M_Sun_kg;    //scaled mass of Earth
  Pos = {1, 0};
  Vel = {0, 2*pi};
  M = ME;
}

Body::Body(vec pos, vec vel, double m)
{
  Pos = pos;      // set internal attribute to given value (e.g. Earth pos)
  Vel = vel;
  M = m;
}

// define member functions for Body
void Body::print()
{
  cout << Pos << '\n' << Vel << '\n' << M << endl;
}

void Body::update(vec newPos, vec newVel)
{
  Pos = newPos;
  Vel = newVel;
  return;
}


double Body::distance(Body otherbody)
{
  double delta_x = Pos(0) - otherbody.Pos(0);
  double delta_y = Pos(1) - otherbody.Pos(1);
//  cout << delta_x << ' ' << delta_y << endl;
  return sqrt(delta_x*delta_x + delta_y*delta_y);
}

double Body::grav_force(Body otherbody)
{
  double d = distance(otherbody);
  double m1 = M;
  double m2 = otherbody.M;
  double F = G*m1*m2/(d*d);

  return F;
}
