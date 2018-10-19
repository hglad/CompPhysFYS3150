#include "Body.h"
#include <armadillo>
// define constructors for Body

Body::Body()
{
  Pos = {-1, 0};
  Vel = {0, -pi};
  M = 1e-6;
}

Body::Body(vec pos, vec vel, double m)
{
  Pos = pos;      // set internal attribute to given value (e.g. Earth pos)
  Vel = vel;
  M = m;
}

// define functions for Body

void Body::print()
{
  cout << Pos << '\n' << Vel << '\n' << M << endl;
}

double Body::distance(Body otherbody)
{
  double delta_x = Pos(0) - otherbody.Pos(0);
  double delta_y = Pos(1) - otherbody.Pos(1);

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
