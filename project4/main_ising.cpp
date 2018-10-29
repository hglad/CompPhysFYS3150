#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <random>

using namespace std;
using namespace arma;

int main(int argc, char const *argv[])
{
  int M = 1000;     // num. of MC-cycles
  int L = 2;  int n = L*L;
  int temp;
  imat S(L,L);

  // Generate L*L matrix with spins -1 or 1
  srand (time(NULL));
  for (int i=0; i < n; i++)
  {
    temp = rand() % 2;    // returns 0 or 1, corresponding to spin -1 and 1
    if (temp == 0)
    {
      S(i) = -1;
    }
    else
    {
      S(i) = 1;
    }
  }

  cout << S << endl;

  

  return 0;

}
