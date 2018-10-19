#ifndef BODY_H
#define BODY_H

#include "project3.h"
#include "math.h"

class Body
{
private:
  double pi = M_PI;
  double G = 4*pi*pi;

public:
  vec Pos;
  vec Vel;
  double M;

  // Constructors: assign given initial values to a body
  Body();        // using generic values

  Body(vec pos, vec vel, double m);    // using specified values

  void print();

  void update(vec newPos, vec newVel);

  double distance(Body otherbody);

  double grav_force(Body otherbody);


};


#endif /* BODY_H */
