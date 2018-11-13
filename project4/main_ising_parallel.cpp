#include "ising.h"
//compile: mpic++ -o3 -o main.x main_ising_parallel.cpp ising.cpp -DARMA_DONT_USE_WRAPPER -lblas -llapack
//run:     mpirun -n 1 ./main.x 20000 2
// Begin main
int main(int argc, char* argv[])
{
  int numMC = atoi(argv[1]);     // num. of MC-cycles
  int L = atoi(argv[2]);
  int n = L*L;
  //int temp_spin;
  double T = atof(argv[2]);

  mat S(L,L);
  double r, energy, magmom;
  double E, E2, M, M2, absM, C_V, chi;
  int x, y, dE;

  vec ValueSums = zeros<vec>(5);              // sum of various parameters
  vec Energy = zeros(numMC); vec Energy2 = zeros(numMC);
  vec Magmom = zeros(numMC); vec Magmom2 = zeros(numMC);
  vec absMagmom = zeros(numMC);
  vec Cv = zeros(numMC);
  vec Chi = zeros(numMC);

  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937_64 gen(10); //Standard mersenne_twister_engine seeded with rd()
  //mt19937_64 gen(rd());
  uniform_real_distribution<double> dist(0.0, 1.0);
  uniform_int_distribution<int> RNGpos(0, L-1);

  // Initial values
  //S.fill(1);              // ordered state
  rand_spins(S);
  init_params(S, energy, magmom);   // initial energy, magnetic momentum

//  myfile << energy << ' ' << magmom << endl;
  int counter;
  map<double, double> w = transitions(T);     // create dictionary

  Energy(0) = energy; Magmom(0) = magmom;
  cout << Energy(0) << endl;
  // Initialize parallellization
  int numprocs, my_rank;
  double initial_temp, final_temp, temp_step;
  double time_start, time_end, total_time;
  double average[5], total[5];
  //double total_energy[numMC];
  vec total_energy = zeros(numMC);

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  time_start = MPI_Wtime();

  // Find proper intervals based on number of processes used
  int no_intervals = numMC/numprocs;
  int mc_start = my_rank*no_intervals + 1;
  int mc_end = (my_rank + 1)*no_intervals;
  if ( (my_rank == numprocs-1) &&( mc_end < numMC) )
  {
    mc_end = numMC;
  }

  // Broadcast variables to allow for parallellization
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //int local_k = my_rank*(numMC/numprocs);   //local start index for processor

  for (int k = mc_start; k < mc_end; k++)
  {
    MC_cycle(S, L, counter, energy, magmom, w, gen);

    Energy(k) = energy;
  //  Energy2(k) = energy*energy;

    Magmom(k) = magmom;
  //  Magmom2(k) = magmom*magmom;

  //  Cv(k) = (Energy2(k) - Energy(k)*Energy(k))/(T*T*n);
  //  Chi(k) = (Magmom2(k) - Magmom(k)*Magmom(k))/(T*n);

  //  myfile << energy << ' ' << magmom << ' ' << Cv(k) << ' ' << Chi(k) << endl;

  //  MPI_Reduce(&Energy[k], &total_energy[k], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    ValueSums(0) += energy; ValueSums(1) += energy*energy;
    ValueSums(2) += magmom; ValueSums(3) += magmom*magmom;
    ValueSums(4) += fabs(magmom);
  }

  // Add all contributions to master node (rank 0)
  for (int i = 0; i < 5; i++)
  {
    MPI_Reduce(&ValueSums(i), &total[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  time_end = MPI_Wtime();
  total_time = time_end - time_start;
  cout << "Time = " <<  total_time  << " on number of processors: "  << numprocs  << endl;
  // Compute expectation values from master node
  if (my_rank == 0)
  {

    write_params(Energy, Magmom);
    //cout << Energy << endl;
    E    = total[0]/numMC;
    E2   = total[1]/numMC;
    M    = total[2]/numMC;
    M2   = total[3]/numMC;
    absM = total[4]/numMC;

    C_V = (E2 - E*E)/(T*T*n);
    chi = (M2 - M*M)/(T*n);
    cout << "---" << endl;
    cout << E << ' ' << M << ' ' << M2 << ' ' << C_V << ' ' << chi << endl;
  }

  MPI_Finalize();

  cout << "Rank:" << ' ' << my_rank << endl;
  //Energy.save("ising_data.txt", arma_ascii);
  return 0;
}
