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

void init_forward(double alpha, double &a, double &c, vec& d);

void init_backward(double alpha, double &a, double &c, vec& d);

void tridiag(double a, double c, vec d, vec& f, vec& u, int n);

#endif /* DIFFUSION_H */
