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
#include <string>

using namespace arma;
using namespace std;

int PBC(int index, int max, int add);

void init_params(mat S, double &E, double &M);

mat init_spins(int L, mt19937_64 &gen, int rand_state);

void write_arrays(vec A, vec B, float T);

void write_means(double E, double absM, double M2, double C_V, double chi, float T);

map<double, double> transitions(double T);

void MC_cycle(mat &S, int L, int& counter, double& energy, double& magmom, map<double, double> w, mt19937_64 &gen);


#endif /* ISING_H */
