#include "ising.h"

int PBC(int index, int max, int add)
{
  return (index + max + add) % max;
}

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
}

void write_params(vec A, vec B)
{
  ofstream myfile;
  myfile.open ("ising_data.txt");
  //cout << A << endl;
  for (int i=0; i < A.n_elem; i++)
  {
    myfile << A(i) << ' ' << B(i) << endl;
  }
  myfile.close();
  return;
}

void rand_spins(mat &S)
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
}

void MC_cycle(mat &S, int L, int& counter, double& energy, double& magmom, map<double, double> w, mt19937_64 &gen)
{     // a single MC cycle
  uniform_int_distribution<int> RNGpos(0, L-1);
  uniform_real_distribution<double> dist(0.0, 1.0);
  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < L; j++)
    {
      int x = RNGpos(gen);
      int y = RNGpos(gen);
      double dE = 2*S(x, y) * (S(x, PBC(y, L, -1)) + S(PBC(x, L, -1), y) + S(x, PBC(y, L, 1)) + S(PBC(x, L, 1), y));
  //    cout << dE << endl;
      // Metropolis algorithm
      // compare w with random number r
      double r = dist(gen);
      if (r <= w.find(dE)->second)    // find corresponding energy to dE
      {
        counter += 1;
        S(x, y)  *= -1;           // flip spin
        energy  += dE;
        magmom  += S(x, y);
      }
    }
  }
  //cout << S(x, y) << endl;
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
