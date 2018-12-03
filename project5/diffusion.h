#ifndef ISING_H
#define ISING_H

#include <armadillo>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <random>

using namespace arma;
using namespace std;

mat init_forward(int n, double alpha);

void tridiag(double a, double c, vec d, vec& f, vec& u, int n);

#endif /* ISING_H */
