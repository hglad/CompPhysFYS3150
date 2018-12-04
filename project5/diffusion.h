#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <armadillo>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <random>

using namespace arma;
using namespace std;

void init_forward(double alpha, double &a, double &c, vec& b);

void init_backward(double alpha, double &a, double &c, vec& b);

void tridiag(double a, double c, vec b, vec& y, vec& u, int nx);

//void write_sol(vec u, string filename);

#endif /* DIFFUSION_H */
