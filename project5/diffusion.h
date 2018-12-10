#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <armadillo>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace arma;
using namespace std;

string init_method(int method, int dim, double dx, double alpha, double &a, double &c, double &b);

void init_backward(double alpha, double &a, double &c, double &b);

void init_crank(double alpha, double &a, double &c, double &b);

void forward(double alpha, vec& u, int nx, double BC1, double BC2);

void backward(double a, double c, double b, double alpha, vec& u, vec y, int nx, double BC1, double BC2);

void crank(double a, double c, double b, double alpha, vec& u, vec y, int nx);

void forward_2D(double alpha, double D, mat& u, int nx, int ny, double BC1, double BC2);

void backward_2D(double a, double c, double b, double alpha, mat& u, mat y, int nx, int ny);

void crank_2D(double a, double c, double b, double alpha, mat& u, mat y, int nx, int ny);

void forward_source(double fac, double fac2, double dt, mat& u, int nx, int ny);

void set_BCs_1D(vec& u, int nx, double BC1, double BC2);

void set_BCs_2D(mat& u, int nx, int ny, double BC1, double BC2);

void tridiag(double a, double c, double b, vec y, vec& u, int nx);

void analytic_1D(int nx, int nt, double L, int saved_steps);

void analytic_2D(int nx, int nt, double L, int saved_steps);

//void write_sol(vec u, string filename);

#endif /* DIFFUSION_H */
