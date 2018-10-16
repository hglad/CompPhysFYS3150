#ifndef BODY_H
#define BODY_H

#include "project3.h"

class Body
{
public:
  vec Pos;
  vec Vel;
  double M;

  // Constructors: assign given initial values to a body
  Body()        // using generic values
  {
    Pos = {-1, 0};
    Vel = {0, -M_PI};
    M = 1e-6;
  }

  Body(vec pos, vec vel, double m)    // using specified values
  {
    Pos = pos;      // set internal attribute to given value (e.g. Earth pos)
    Vel = vel;
    M = m;

  }


};


/*
int main()
{
  double pi = M_PI;
  vec s_E = {1, 0};
  vec v_E = {0, 2*pi};
  double m_E = 6*pow(10, 24) / 2*pow(10, 30);

  Body planet(s_E, v_E, m_E);    // instance of 'body'
  return 0;

}
*/

#endif /* BODY_H */
