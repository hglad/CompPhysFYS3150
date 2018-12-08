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

void forward(double alpha, vec& u, int nx);

void backward(double a, double c, double b, double alpha, vec& u, vec y, int nx);

void crank(double a, double c, double b, double alpha, vec& u, vec y, int nx);

void forward_2D(double alpha, mat& u, int nx, int ny);

void backward_2D(double a, double c, double b, double alpha, mat& u, mat y, int nx, int ny);

void crank_2D(double a, double c, double b, double alpha, mat& u, mat y, int nx, int ny);

void forward_source(double fac, double fac2, double dt, mat& u, int nx, int ny);

void set_BCs_2D(mat& u);

void set_BCs_source(mat& u);

void tridiag(double a, double c, double b, vec y, vec& u, int nx);

void analytic(int nx, int nt);

//void write_sol(vec u, string filename);

#endif /* DIFFUSION_H */
