#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

vec arma_diag(mat A);
double max_offdiag(mat A, int *k, int *l, double& max_elem);
void Jacobi_rot(mat& A, double& max_elem);


#endif /* JACOBI_H */
