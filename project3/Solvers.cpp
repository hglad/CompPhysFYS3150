#include "Solvers.h"
#include "Body.h"
#include <vector>

// if no parameters are specified, calculate orbit of Earth with
// only the Sun as a source of attraction
Solvers::Solvers()
{
  dt = 0.01;
  Body target_Body();     // use default constructor in Body class
  // vector with only the sun as a Body object
  vec zero = {0,0};
  Body sun(zero, zero, 1);
  Bodies = { sun };
  n_ext_bodies = 1;

}

Solvers::Solvers(Body target_body, vector<Body> bodies, double DT)
{
  dt = DT;
  target_Body = target_body;
  Bodies = bodies;
  n_ext_bodies = bodies.size();
}

// define functions for Solvers class

void Solvers::Verlet(vec& pos, vec& vel, vec& old_Acc)
{
  // Define old values of position, velocity, acceleration - needed for Verlet
  old_Pos = pos;
  old_Vel = vel;
  pos(0) = old_Pos(0) + dt*old_Vel(0) + old_Acc(0)*dt*dt/2;
  pos(1) = old_Pos(1) + dt*old_Vel(1) + old_Acc(1)*dt*dt/2;

  // Loop through bodies in the system that accelerate the target body
  for (int j=0; j < n_ext_bodies; j++)
  {
    temp_d = target_Body.distance(Bodies[j]);
    cout << temp_d << endl;
    temp_a = target_Body.grav_force(Bodies[j])/target_Body.M;
    delta_x = target_Body.Pos(0) - Bodies[j].Pos(0);
    delta_y = target_Body.Pos(1) - Bodies[j].Pos(1);
    acc(0) += -temp_a * delta_x/temp_d;     // x-component of acceleration
    acc(1) += -temp_a * delta_y/temp_d;     // y-component of acceleration
  }

  vel(0) = old_Vel(0) + dt/2*( acc(0) + old_Acc(0) );
  vel(1) = old_Vel(1) + dt/2*( acc(1) + old_Acc(1) );

  old_Acc = acc;

  return;
}
