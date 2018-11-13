#ifndef ISING_H
#define ISING_H

#include <armadillo>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <random>
#include <map>
#include "mpi.h"

using namespace arma;
using namespace std;

int PBC(int index, int max, int add);

void init_params(mat S, double &E, double &M);

void write_params(vec A, vec B);

void rand_spins(mat &S);

map<double, double> transitions(double T);

void MC_cycle(mat &S, int L, int& counter, double& energy, double& magmom, map<double, double> w, mt19937_64 &gen);


#endif /* ISING_H */
