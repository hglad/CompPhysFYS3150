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

void write_arrays(vec A, vec B, int numMC, int L, float T)
{
  string temp = to_string(T);
  temp = temp.substr(0,4);
  string dim = to_string(L);
  string MC = to_string(numMC);
  ofstream myfile;
  myfile.open ("results/arrays_L=" + dim + "T=" + temp + "MC=" + MC + ".txt");

  for (int i=0; i < numMC; i++)
  {
    myfile << A(i) << ' ' << B(i) << endl;
  }

  return;
}

void write_means(vec E, vec absM, vec M2, vec C_V, vec chi, int *counter, int numMC, int L, vec T_vec)
{
  string dim = to_string(L);
  string MC = to_string(numMC);

  int n_temps = T_vec.n_elem;
  double T_start = T_vec(0);
  double T_final = T_vec(n_temps-1);

  ofstream myfile;
  string temp_start = to_string(T_start);
  temp_start = temp_start.substr(0,4);

  string temp_final = to_string(T_final);
  temp_final = temp_final.substr(0,4);
  myfile.open ("results/means_L=" + dim + "T=" + temp_start + "-" + temp_final + "MC=" + MC + ".txt");

  for (int i=0; i < n_temps; i++)
  {
    myfile << E(i) << ' ' << absM(i) << ' ' << M2(i) << ' ' << C_V(i) << ' ' << chi(i) << ' ' << counter[i] << ' ' << T_vec(i) << ' ' << MC << endl;
  }
  myfile.close();
  return;
}

// Generate L*L matrix with spins -1 or 1
mat init_spins(int L, mt19937_64 &gen, int rand_state)
{
  mat S(L,L);
  if (rand_state == 1)
  {
    uniform_int_distribution<int> RNGspin(0, 1);
    int n = L*L;
    for (int i=0; i < n; i++)
    {
      S(i) = RNGspin(gen);  // returns 0 or 1, corresponding to spin -1 and 1
      if (S(i) == 0)
      {
        S(i) = -1;
      }
    }
  }
  else
  {
    S.fill(1);
  }
return S;
}

void MC_cycle(mat &S, int L, int& counter, double& energy, double& magmom, map<double, double> w, mt19937_64 &gen)
{     // a single MC cycle
  int n = L*L;
  uniform_int_distribution<int> RNGpos(0, L-1);
  uniform_real_distribution<double> dist(0.0, 1.0);
  for (int i = 0; i < n; i++)
  {
      int x = RNGpos(gen);
      int y = RNGpos(gen);

      double dE = 2*S(x, y) * (S(x, PBC(y, L, -1)) + S(PBC(x, L, -1), y) + S(x, PBC(y, L, 1)) + S(PBC(x, L, 1), y));
      // Metropolis algorithm:
      // compare w with random number using dist(gen)
      if (dist(gen) <= w.find(dE)->second)   // find corresponding energy to dE
      {
        counter += 1;            // counts number of states accepted
        S(x, y) *= -1;           // flip spin
        energy  += dE;
        magmom  += 2*S(x, y);    // difference in mag. after flip
      }
  }
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
// Function to reset sums needed for expectation values inbetween calculations
void reset_sums(vec &ValueSums, vec &total)
{
  ValueSums(0) = 0; ValueSums(1) = 0;
  ValueSums(2) = 0; ValueSums(3) = 0;
  ValueSums(4) = 0;

  total(0) = 0; total(1) = 0;
  total(2) = 0; total(3) = 0;
  total(4) = 0;
  return;
}
