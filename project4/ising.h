#ifndef ISING_H
#define ISING_H

#include <armadillo>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <random>
#include <map>
#include <string>

using namespace arma;
using namespace std;

int PBC(int index, int max, int add);

void init_params(mat S, double &E, double &M);

mat init_spins(int L, mt19937_64 &gen, int rand_state);

void write_arrays(vec A, vec B, int numMC, int L, float T);

void write_means(vec E, vec absM, vec M2, vec C_V, vec chi, int *counter, int numMC, int L, vec T_vec);

map<double, double> transitions(double T);

void MC_cycle(mat &S, int L, int& counter, double& energy, double& magmom, map<double, double> w, mt19937_64 &gen);

void reset_sums(vec &ValueSums, vec &total);

#endif /* ISING_H */
