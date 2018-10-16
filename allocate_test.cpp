#define _USE_MATH_DEFINES

#include <iostream>
//#include <armadillo>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
  int n=100;
  double **A;
  A = new double*[n];

  for (int i=0; i < n; i++)
  {
    A[i] = new double[n];
  }

  cout << A[0][0][0] << endl;

  return 0;
}
