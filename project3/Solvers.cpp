#include "Solvers.h"
#include "Body.h"
#include <vector>
#include <fstream>

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
  // Define old values of position - needed for Verlet
  old_Pos = pos;
  old_Vel = vel;
  acc = {0,0};    // reset acceleration

  pos = old_Pos + dt*old_Vel + old_Acc*dt*dt/2;
  // Loop through bodies in the system that accelerate the target body
  for (int j=0; j < n_ext_bodies; j++)
  {
    temp_d = target_Body.distance(Bodies[j]);
    cout << temp_d << endl;
    temp_a = target_Body.grav_force(Bodies[j])/target_Body.M; // acceleration magnitude

    Delta_x = target_Body.Pos(0) - Bodies[j].Pos(0); // x-coordinate difference
    Delta_y = target_Body.Pos(1) - Bodies[j].Pos(1); // y-coordinate difference

    acc(0) += -temp_a * Delta_x/temp_d;     // x-component of acceleration
    acc(1) += -temp_a * Delta_y/temp_d;     // y-component of acceleration
  }

  vel = old_Vel + dt/2*( acc + old_Acc );
  old_Acc = acc;

  target_Body.update(pos, vel);

  return;
}

void Solvers::Euler(vec& pos, vec& vel)
{
  acc = {0,0};    // reset acceleration
  for (int j=0; j < n_ext_bodies; j++)
  {
    temp_d = target_Body.distance(Bodies[j]);
    cout << temp_d << endl;
    temp_a = target_Body.grav_force(Bodies[j])/target_Body.M; // acceleration magnitude

    Delta_x = target_Body.Pos(0) - Bodies[j].Pos(0); // x-coordinate difference
    Delta_y = target_Body.Pos(1) - Bodies[j].Pos(1); // y-coordinate difference

    acc(0) += -temp_a * Delta_x/temp_d;     // x-component of acceleration
    acc(1) += -temp_a * Delta_y/temp_d;     // y-component of acceleration
  }

  vel = vel + acc*dt;
  pos = pos + vel*dt;

  target_Body.update(pos, vel);

  return;
}
