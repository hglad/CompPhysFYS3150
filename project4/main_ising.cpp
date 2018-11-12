#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <random>
#include <fstream>
#include <map>
#include <mpi.h>

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
  int numMC = 100000;     // num. of MC-cycles
  int L = 2;  int n = L*L;
  //int temp_spin;
  double T = 1;

  imat S(L,L);
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

  ofstream myfile;
  myfile.open ("ising_data.txt");
  myfile << energy << ' ' << magmom << endl;

  int counter;
  map<double, double> w = transitions(T);     // create dictionary
  vec ExpectVals = zeros<vec>(5);

  // Parallellization
  int numprocs, my_rank;
  double initial_temp, final_temp, temp_step;
  double time_start, time_end, total_time;
  double average[5], total[5];

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  time_start = MPI_Wtime();

  // Find number of MC cycles based on number of processes used
  int no_intervalls = numMC/numprocs;
  int myloop_begin = my_rank*no_intervalls + 1;
  int myloop_end = (my_rank+1)*no_intervalls;
  if ( (my_rank == numprocs-1) &&( myloop_end < numMC) ) myloop_end = numMC;

  // Broadcast variables to allow for parallellization
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  int local_n = n/numprocs;

  //for (int k=1; k <= numMC; k++)
  for (int k = myloop_begin; k <= myloop_end; k++)
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        x = RNGpos(gen);
        y = RNGpos(gen);
        dE = 2*S(x, y) * (S(x, PBC(y, L, -1)) + S(PBC(x, L, -1), y) + S(x, PBC(y, L, 1)) + S(PBC(x, L, 1), y));

        // Metropolis algorithm
        // compare w with random number r
        r = dist(gen);
        if (r <= w.find(dE)->second)    // find corresponding energy to dE
        {
          counter += 1;
          S(x, y) *= -1;           // flip spin
          energy  += dE;
          magmom  += S(x, y);
          }
      }
    }
    /*
    if (my_rank == 0)
    {
      myfile << energy << ' ' << fabs(magmom) << endl;
    }
    */
    ExpectVals(0) += energy; ExpectVals(1) += energy*energy;
    ExpectVals(2) += magmom; ExpectVals(3) += magmom*magmom;
    ExpectVals(4) += fabs(magmom);
  }

  // Compute expectation values
  for (int i = 0; i < 5; i++)
  {
    MPI_Reduce(&ExpectVals[i], &total[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (my_rank == 0)
  {
    for (int i = 0; i < 5; i++)
    {
      cout << total[i]/(numMC) << endl;
    }

  }
  /*
  E  = ExpectVals(0)/numMC;
  E2 = ExpectVals(1)/numMC;
  M  = ExpectVals(2)/numMC;
  M2 = ExpectVals(3)/numMC;
  */
  C_V = (E2 - E*E)/(pow(T,2)*pow(L,2));
  chi = (M2 - M*M)/(T*pow(L,2));

//  cout << ExpectVals/numMC << endl;
  // cout << "---" << endl;
//  cout << E << ' ' << M << ' ' << M2 << ' ' << C_V << ' ' << chi << endl;

  time_end = MPI_Wtime();
  total_time = time_end - time_start;
  if ( my_rank == 0)
  {
    cout << "Time = " <<  total_time  << " on number of processors: "  << numprocs  << endl;
  }

  myfile.close();
  MPI_Finalize ();

  return 0;

}
/*
  void flip(mat& S)        // function to flip a random spin
  {
    int n = S.n_elem;
    int rand_pos = rand() % n;
    S(rand_pos) *= -1;      // changing sign --> flipping spin
    return;
  }
*/
