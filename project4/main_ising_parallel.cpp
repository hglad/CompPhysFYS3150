#include "ising.h"
//compile: mpic++ -o3 -o main.x main_ising_parallel.cpp ising.cpp -DARMA_DONT_USE_WRAPPER -lblas -llapack
//run:     mpirun -n 1 ./main.x 20000 2 1 2.4 1
//arguments: numMC L T_start T_end rand_state cut_off_num

int main(int argc, char* argv[])
{
  // ------ Initialize variables ------
  int numMC = atoi(argv[1]);     // num. of MC-cycles
  int L = atoi(argv[2]);
  double T_start = atof(argv[3]);
  double T_final = atof(argv[4]);
  int rand_state = atoi(argv[5]); // rand_state=1 gives a random initial configuration
  int cut_off_num = atoi(argv[6]); // number of initial MC-cycles to cut off
  int n = L*L;

  double energy, magmom;

  bool save_arrays = false;
  bool save_means = true;

  //cout << vec_len <<endl;

  vec ValueSums = zeros<vec>(5);              // sum of various parameters
  //vec Energy2 = zeros(numMC);
  //vec Magmom2 = zeros(numMC);
  //vec absMagmom = zeros(numMC);
  //vec Cv = zeros(numMC);
  //vec Chi = zeros(numMC);

  // ------ Initialize random number generator ------
  random_device rd;  //Will be used to obtain a seed for the random number engine
  //mt19937_64 gen(10); //Standard mersenne_twister_engine seeded with rd()
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

  // ------ Set up initial spins and initial energy, magnetization ------
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
//  ivec counter = zeros<ivec>(n_temps);
  int *counter = new int[n_temps];

  int i = 0;    // counter for number of temperatures computed
  int start_sum = mc_start + cut_off_num;
//  int counter__;
  //int local_k = my_rank*(numMC/numprocs);   //local start index for processor
  for (double T = T_start; T <= T_final*1.0001; T+=T_step)
  {
    counter[i] = 0;
    time_start = MPI_Wtime();
    map<double, double> w = transitions(T);     // create dictionary
    //counter = 0;              // accepted states
    reset_sums(ValueSums, total);   // reset sums for expectation values
    for (int k = mc_start; k < mc_end; k++)
    {
      MC_cycle(S, L, counter[i], energy, magmom, w, gen);

      if (my_rank == 0)
      {
        Energy(k) = energy;
        Magmom(k) = fabs(magmom);
      }

      if (k > start_sum)     // do not use first 3000 cycles for every process
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
      //cut_off_mc = numMC - 3000;
      total /= (numMC - cut_off_num*numprocs);
      // multiply cut off with number of processors since every process cuts off
      // the same amount of cycles

      E(i)    = total(0);
      E2(i)   = total(1);
      M(i)    = total(2);
      M2(i)   = total(3);
      absM(i) = total(4);

      C_V(i) = (E2(i) - E(i)*E(i))/(T*T);
      chi(i) = (M2(i) - absM(i)*absM(i))/(T);

      cout << "Results: T = " << T << endl;
    //  cout << E(i)/n << ' ' << absM(i)/n << ' ' << M2(i)/n << ' ' << C_V(i)/n << ' ' << chi(i)/n << endl;
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

  //Energy.save("ising_data.txt", arma_ascii);
  return 0;
}

// send/receive routine (unused)
/*
if (my_rank != 0)
{
  MPI_Send(Energy, numMC, MPI_DOUBLE, 0, 50, MPI_COMM_WORLD);
}

if (my_rank == 0)
{
  double *energies = new double[numMC];
  string temp = to_string(T);
  MPI_Status status;
  ofstream myfile;
  // open and close to clear previous contents of file
  myfile.open ("ising_arrays_" + temp + ".txt");
  myfile.close();

  myfile.open ("ising_arrays_" + temp + ".txt", ios_base::app);

  for (int i=0; i < numMC; i++)
  {
      myfile << energies[i] << ' ' << Magmom(i) << endl;
  }

  for (int rank=1; rank < numprocs; rank++)
  {
    MPI_Recv(energies, numMC, MPI_DOUBLE, rank, 50, MPI_COMM_WORLD, &status);

    for (int i=0; i < numMC; i++)
    {
      myfile << energies[i] << ' ' << Magmom(i) << endl;
    }

    //write_arrays(energies, Magmom, numMC, T);
  }
  myfile.close();

}
*/
