/* -------- Documentation --------

compile: mpic++ -o3 -o main.x main_ising_parallel.cpp ising.cpp -DARMA_DONT_USE_WRAPPER -lblas -llapack

command line arguments in order:
numMC - total number of Monte Carlo cycles to perform (independent of processors)
L - size of lattice dimensions
T_start - first temperature to calculate expectation values for
T_end - last temperature (temperature step size is 0.01)
rand_state - controls if initial configuration is random or ordered, set to 1 for random
cut_off_num - number of initial MC-cycles to cut off per processor

Functions are located in 'ising.cpp'. By default, running the program does not write
individual energies and magnetic momentum to file. This is to save computation time
but can be enabled by setting the boolean variable 'save_arrays' to 'true'.

IMPORTANT: make sure that a 'results' folder already exists in the same location as
this file, as the values produced are written to a folder with this name. This is
to avoid clutter in the project folder, as a lot of data is produced.
*/
#include "ising.h"
#include "mpi.h"
int main(int argc, char* argv[])
{

  if (argc < 7)
	{
		cout << "Input error: too few arguments. See file documentation." << endl;
		exit(1);
	}
	if (argc > 7)
	{
		cout << "Input error: too many arguments. See file documentation." << endl;
		exit(1);
	}

  // Initialize variables
  int numMC = atoi(argv[1]);
  int L = atoi(argv[2]);
  double T_start = atof(argv[3]);
  double T_final = atof(argv[4]);
  int rand_state = atoi(argv[5]);
  int cut_off_num = atoi(argv[6]);
  int n = L*L;

  double energy, magmom;

  bool save_arrays = false;
  bool save_means = true;

  vec ValueSums = zeros<vec>(5);              // sum of various parameters

  // Initialize random number generators
  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937_64 gen(rd());
  uniform_real_distribution<double> dist(0.0, 1.0);
  uniform_int_distribution<int> RNGpos(0, L-1);

  // Initialize parallellization
  int numprocs, my_rank;
  double time_start, time_end, total_time;
  vec total = zeros<vec>(5);

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  // Find proper intervals based on number of processes used
  int no_intervals = numMC/numprocs;

  vec Energy = zeros(no_intervals);
  vec Magmom = zeros(no_intervals);
  int mc_start = my_rank*no_intervals + 1;
  int mc_end = (my_rank + 1)*no_intervals;
  if ( (my_rank == numprocs-1) &&( mc_end < numMC) )
  {
    mc_end = numMC;
  }

  // Set up initial spins and initial energy, magnetization
  mat S = init_spins(L, gen, rand_state);
  init_params(S, energy, magmom);   // initial energy, magnetic momentum
  Energy(0) = energy; Magmom(0) = fabs(magmom);

  // Broadcast variables
  double T_step = 0.01;
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&T_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&T_final, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&T_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Create vectors to hold expectation values for different temperatures
  int n_temps = ceil((T_final - T_start)/T_step) + 1;

  vec E = zeros(n_temps); vec E2 = zeros(n_temps);
  vec M = zeros(n_temps); vec M2 = zeros(n_temps);
  vec absM = zeros(n_temps); vec C_V = zeros(n_temps);
  vec chi = zeros(n_temps);
  int *counter = new int[n_temps];

  int i = 0;    // counter for number of temperatures computed
  int start_sum = mc_start + cut_off_num;
  for (double T = T_start; T <= T_final*1.0001; T+=T_step)
  {
    counter[i] = 0;
    time_start = MPI_Wtime();
    map<double, double> w = transitions(T); // create dictionary
    reset_sums(ValueSums, total);           // reset sums for expectation values
    for (int k = mc_start; k < mc_end; k++)
    {
      MC_cycle(S, L, counter[i], energy, magmom, w, gen);

      if (my_rank == 0)
      {
        Energy(k) = energy;
        Magmom(k) = fabs(magmom);
      }

      if (k > start_sum)     // do not add first MC-cycles to sum
      {
        ValueSums(0) += energy; ValueSums(1) += energy*energy;
        ValueSums(2) += magmom; ValueSums(3) += magmom*magmom;
        ValueSums(4) += fabs(magmom);
      }
    }

    // Add all contributions to master node (rank 0)
    for (int j = 0; j < 5; j++)
    {
      MPI_Reduce(&ValueSums(j), &total(j), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Compute expectation values from master node, write values
    if (my_rank == 0)
    {
      time_end = MPI_Wtime();
      total_time = time_end - time_start;
      cout << "Time = " <<  total_time << " on number of processors: " << numprocs << endl;
      total /= (numMC - cut_off_num*numprocs);
      // multiply cut off with number of processors since every process cuts off
      // the same amount of cycles

      E(i)    = total(0);
      E2(i)   = total(1);
      M(i)    = total(2);
      M2(i)   = total(3);
      absM(i) = total(4);

      C_V(i) = (E2(i) - E(i)*E(i))/(T*T);
      chi(i) = (M2(i) - M(i)*M(i))/(T);

      cout << "Results: T = " << T << endl;
      printf ("E = %1.4f absM = %1.4f M**2 = %1.4f C_V = %1.4f chi = %1.4f\n", E(i), absM(i), M2(i), C_V(i), chi(i));

      if (save_arrays == true)
      {
        write_arrays(Energy, Magmom, no_intervals, L, T);
      }
      i += 1;
      // only write arrays from one process
    }

  }
  // save expectation values from all different temperatures
  if ((save_means == true) && (my_rank == 0))
  {
    vec T_vec = linspace(T_start, T_final, n_temps);
    write_means(E/n, absM/n, M2/n, C_V/n, chi/n, counter, numMC, L, T_vec);
  }

  MPI_Finalize();

  return 0;
}
