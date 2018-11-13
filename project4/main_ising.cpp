#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <random>
#include <fstream>
#include <map>

using namespace std;
using namespace arma;

inline double energy_calc(mat S) {return -(2*S(0)*S(1) + 2*S(0)*S(2) + 2*S(1)*S(3) + 2*S(2)*S(3));}

inline int PBC(int index, int max, int add) {return (index + max + add) % max ;}

void init_params(mat S, double &E, double &M)
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
  return;
}

void write_params(vec A, vec B)
{
  ofstream myfile;
  myfile.open ("ising_data.txt");

  for (int i=0; i < A.n_elem; i++)
  {
    myfile << A(i) << ' ' << B(i) << endl;
  }
  myfile.close();
  return;
}

void rand_spins(imat &S)
{
  // Generate L*L matrix with spins -1 or 1
  srand (time(NULL));
  int n = S.n_elem;
  for (int i=0; i < n; i++)
  {
    S(i) = rand() % 2;  // returns 0 or 1, corresponding to spin -1 and 1
    if (S(i) == 0)
    {
      S(i) = -1;
    }
  }
  return;
}

void MC_cycle(int& test, mat &S, int L, int& counter, double& energy, double& magmom, map<double, double> w, int x, int y, double r)
{     // a single MC cycle
  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < L; j++)
    {
      double dE = 2*S(x, y) * (S(x, PBC(y, L, -1)) + S(PBC(x, L, -1), y) + S(x, PBC(y, L, 1)) + S(PBC(x, L, 1), y));
  //    cout << dE << endl;
      // Metropolis algorithm
      // compare w with random number r
      if (r <= w.find(dE)->second)    // find corresponding energy to dE
      {
        counter += 1;
        S(x, y)  *= -1;           // flip spin
        energy  += dE;
        magmom  += S(x, y);
      }
    }
  }
  cout << S(x, y) << endl;
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


// Begin main
int main(int argc, char* argv[])
{
  int numMC = atoi(argv[1]);     // num. of MC-cycles
  int L = atoi(argv[2]);
  int n = L*L;
  //int temp_spin;
  double T = atof(argv[3]);

  mat S(L,L);
  double r, energy, magmom;
  double E, E2, M, M2, absM, C_V, chi;
  int x, y, dE;

  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937_64 gen(10); //Standard mersenne_twister_engine seeded with rd()
  //mt19937_64 gen(rd());
  uniform_real_distribution<double> dist(0.0, 1.0);
  uniform_int_distribution<int> RNGpos(0, L-1);

  // Initial values
  S.fill(1);              // ordered state
  init_params(S, energy, magmom);

//  myfile << energy << ' ' << magmom << endl;

  int counter = 0;
  map<double, double> w = transitions(T);     // create dictionary
  vec ValueSums = zeros<vec>(5);              // sum of various parameters
  vec TotalSums = zeros<vec>(5);

  vec Energy = zeros(numMC); vec Energy2 = zeros(numMC);
  vec Magmom = zeros(numMC); vec Magmom2 = zeros(numMC);
  vec absMagmom = zeros(numMC);
  Energy(0) = energy; Magmom(0) = magmom;

  ofstream myfile;
  myfile.open ("ising_data.txt");
  myfile << energy << ' ' << magmom << endl;

  for (int k = 1; k < numMC; k++)
  {
    x = RNGpos(gen);
    y = RNGpos(gen);
    r = dist(gen);
    MC_cycle(S, L, counter, energy, magmom, w, x, y, r);

    myfile << energy << ' ' << magmom << endl;
    //cout << energy << endl;
    //cout << S << endl;

    ValueSums(0) += energy; ValueSums(1) += energy*energy;
    ValueSums(2) += magmom; ValueSums(3) += magmom*magmom;
    ValueSums(4) += fabs(magmom);
  } //end for

  cout << counter << endl;
  // Compute expectation values (mean)

    E    = ValueSums(0)/numMC;
    E2   = ValueSums(1)/numMC;
    M    = ValueSums(2)/numMC;
    M2   = ValueSums(3)/numMC;
    absM = ValueSums(4)/numMC;

    C_V = (E2 - E*E)/(T*T*n);
    chi = (M2 - M*M)/(T*n);

    cout << "---" << endl;
    cout << E << ' ' << M << ' ' << M2 << ' ' << C_V << ' ' << chi << endl;

  //write_params(Energy, Magmom);
  //Energy.save("ising_data.txt", arma_ascii);
  myfile.close();
  return 0;
}
