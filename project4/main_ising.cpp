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

inline int PBC(int index, int max, int add) {return (index + max + add) % max ;}

void init_params(imat S, double &E, double &M)
{
  int L = S.n_cols;
  for (int x=0; x < L; x++)
  {
    for (int y=0; y < L; y++)
    {
      M += S(x, y);
      E -= S(x, y)*( S(PBC(x, L, -1), y) + S(x, PBC(y, L, -1)));
    }
  }
}

void newPos(int &x, int &y, int L)
{
  x = rand() % L;
  y = rand() % L;
  return;
}

map<double, double> transitions(double T)
{
  map<double, double> possible_E; // similar to python dictionary

  for (int dE = -8; dE <= 8; dE += 4)
  {
    // dictionary element where dE is the key to the corresponding energy
    possible_E.insert(pair<double, double>(dE, exp(-1/T*(dE))));
  }

  return possible_E;
}

int main(int argc, char const *argv[])
{
  int numMC = 100000;     // num. of MC-cycles
  int L = 2;  int n = L*L;
  //int temp_spin;
  double T = 1;

  imat S(L,L);
  double dE, r, M;
  int x, y;

  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937_64 gen(10); //Standard mersenne_twister_engine seeded with rd()
  //mt19937_64 gen(rd());
  uniform_real_distribution<double> dist(0, 1);
  uniform_int_distribution<int> RNGpos(0, L-1);

  // Generate L*L matrix with spins -1 or 1
  srand (time(NULL));
  for (int i=0; i < n; i++)
  {
    S(i) = rand() % 2;  // returns 0 or 1, corresponding to spin -1 and 1
    if (S(i) == 0)
    {
      S(i) = -1;
    }
  }
  S.fill(1);

  // Initial energy
  double E_0 = 0;
  double magmom = 0;
  init_params(S, E_0, magmom);
  cout << E_0 << endl;
  double energy = E_0;

  ofstream myfile;
  myfile.open ("ising_data.txt");
  myfile << E_0 << ' ' << magmom << endl;

  int counter;
  map<double, double> w = transitions(T);

  for (int k=1; k < numMC; k++)
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
      //  newPos(x, y, L);   // new random position in lattice
        x = RNGpos(gen);
        y = RNGpos(gen);
        dE = 2*S(x, y) * (S(x, PBC(y, L, -1)) + S(PBC(x, L, -1), y) + S(x, PBC(y, L, 1)) + S(PBC(x, L, 1), y));

        // Metropolis algorithm
        // compare w with random number r
        r = dist(gen);
        if (r <= w.find(dE)->second)
        {
          counter += 1;
          S(x, y) *= -1.;           // flip spin
          energy += dE;
          cout << energy << endl;
          magmom += S(x, y);
          }
      }
    }

    myfile << energy << ' ' << fabs(magmom) << endl;
  }
  double E, E2, M2, absM, C_V, chi;
  E = energy/numMC;
  E2 = energy*energy/numMC;
  M = magmom/numMC;
  absM = abs(magmom)/numMC;
  M2 = magmom*magmom/numMC;

  C_V = (E2 - E*E)/(pow(T,2)*pow(L,2));
  chi = (M2 - M*M)/(T*pow(L,2));

  //cout << counter << endl;
  cout << "---" << endl;
  cout << E << ' ' << M << ' ' << M2 << ' ' << C_V << ' ' << chi << endl;

  myfile.close();

  return 0;

}

/*
  void flip(imat& S)        // function to flip a random spin
  {
    int n = S.n_elem;
    int rand_pos = rand() % n;
    S(rand_pos) *= -1;      // changing sign --> flipping spin
    return;
  }
*/
