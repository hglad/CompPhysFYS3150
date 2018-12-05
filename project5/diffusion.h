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

string init_method(int method, int dim, double dx, double alpha, double &a, double &c, vec &b);

void init_backward(double alpha, double &a, double &c, vec& b);

void init_crank(double alpha, double &a, double &c, vec& b);

void forward(double alpha, vec& u, int nx);

void backward(double a, double c, vec b, double alpha, vec& u, vec y, int nx);

void crank(double a, double c, vec b, double alpha, vec& u, vec y, int nx);

void forward_2D(double alpha, mat& u, int nx, int ny);

void backward_2D(double a, double c, vec b, double alpha, mat& u, mat y, int nx, int ny);

void crank_2D(double a, double c, vec b, double alpha, mat& u, mat y, int nx, int ny);

void set_BCs_2D(mat& u);

void tridiag(double a, double c, vec b, vec y, vec& u, int nx);

//void write_sol(vec u, string filename);

#endif /* DIFFUSION_H */