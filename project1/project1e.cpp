#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
using namespace arma;

int main(int argc, char* argv[])
{
  int n=atoi(argv[1]);
//  double *A = new double[n*n];
//  mat arma_A(A, n, n, true, true);

  mat B = zeros(n,n);
  for (int i=0; i < n; i++)
  {
    B[i, i+1]  = -1;
    B[i+1, i]  = -1;
//    B[i+1, i]  = -1;
//    B[i, i+1]  = -1;
    B[i, i]    =  2;
  }

  for (int i=0; i < n; i++)
  {
    std::cout << B[i+1,i] << endl;
  }
  //lu(a);
  return 0;
}
