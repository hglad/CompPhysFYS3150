#include "project3.h"

// Calculate gravitational force on a body, given all bodies in the system
// take masses and coordinate vectors of planets, index j representing target planet
// s(planet, x), s(planet, y)
void grav_force(vec M, mat s, mat& a, int j)
{
  int np = M.n_elem;
  double G = 4*M_PI*M_PI;
  for (int i=0; i < np; i++)
  {
    // Difference between (x,y) coordinates of the two planets
    double delta_x = s(j,0)-s(i,0);
    double delta_y = s(j,1)-s(i,1);
    // Distance planet - planet
    double r = sqrt( (delta_x*delta_x) + (delta_y*delta_y) );
    double F = G*M(i)*M(j)/(r*r);

    a(j,0) -= F/M(j) * delta_x/r;    //cos(theta) = dx/r
    a(j,1) -= F/M(j) * delta_y/r;
  }
  return;
}
