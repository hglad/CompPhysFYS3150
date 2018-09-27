#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

vec arma_diag(mat A);
double max_offdiag(mat A, int *k, int *l, double& max_elem);
void Jacobi_rot(mat& A, double& max_elem);
mat generate_mat(int m, vec d, double a1, double a2);

#endif /* JACOBI_H */
