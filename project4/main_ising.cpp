#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <random>
#include <fstream>
#include <map>

using namespace std;
using namespace arma;

inline double energy_calc(imat S) {return -(2*S(0)*S(1) + 2*S(0)*S(2) + 2*S(1)*S(3) + 2*S(2)*S(3));}

void newPos(int *x, int *y, int L)
{
  *x = rand() % L;
  *y = rand() % L;
  return;
}

void flip(imat& S)        // function to flip a random spin
{
  int n = S.n_elem;
  int rand_pos = rand() % n;
  S(rand_pos) *= -1;      // changing sign --> flipping spin
  return;
}

map<double, double> transitions(double T)
{
  map<double, double> acceptAmp; // similar to python dictionary

  for (int dE = -8; dE <= 8; dE += 4)
  {
    // dictionary element where dE is the key to the corresponding energy
    acceptAmp.insert(pair<double, double>(dE, exp(-1/T*(dE))));
  }

  return acceptAmp;
}

int main(int argc, char const *argv[])
{
  int numMC = 10000;     // num. of MC-cycles
  int L = 2;  int n = L*L;
  int temp_spin;
  double T = 1;

  imat S(L,L);
  vec E = zeros(numMC);
  vec Cvs = zeros(numMC);

  double dE, Cv, M, new_E, w, r;
  int x, y, rand_pos;

  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<> dist(0, 1);

  // Generate L*L matrix with spins -1 or 1
  srand (time(NULL));
  for (int i=0; i < n; i++)
  {
    S(i) = rand() % 2;  // returns 0 or 1, corresponding to spin -1 and 1
    if (temp_spin == 0)
    {
      S(i) = -1;
    }
  }

  cout << S << endl;

  // Initial energy
  double E_0 = energy_calc(S);
  E(0) = E_0;

  ofstream myfile;
  myfile.open ("ising_data.txt");

  for (int i=1; i < numMC; i++)
  {
    flip(S);
    E(i) = energy_calc(S);
    dE = E(i) - E(i-1);

    // Metropolis algorithm
    if (dE > 0)
    {
      // compare w with random number r
      w = exp(-dE/T);
      r = dist(gen);

      // keep new configuration if true
      if (r <= w)
      {
        E(i) = E(i);
      }
      else
      {
        E(i) = E(i-1);
      }

    }

    myfile << E(i) << ' ' << E(i) << endl;
  }

  cout << S << endl;

  myfile.close();

  return 0;

}

/*
flip(S);

new_E = energy(S);
dE = new_E - E;

if (dE <= 0)
{
  E = new_E;
}

else
{
  // compare w with random number r
  w = exp(-dE/T);
  r = dis(gen);
  // keep new configuration if true
  if (r <= w)
  {
    E = new_E;
  }

}
*/
