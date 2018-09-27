#include "project2.h"

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]);
  double rho_max = atoi(argv[2]);
  double rho_0 = 0;
  double h = rho_max/(n);    // rho_min = 0

  //vec rho = linspace(rho_0 + h, rho_max - h, n);
  vec rho = linspace(1, n, n)*h;
  double omega_r = atof(argv[3]);
  vec d = 2./(h*h) + omega_r*omega_r*(rho%rho) + 1./rho;

  double a = -1./(h*h);
  // Generate tridiagonal matrix with added potential
  mat A = generate_mat(n, d, a, a);

  //cout << A << endl;
  double eps = 1e-8;
  double max_elem = 2*eps; // initialize max variable with arbitrary value > eps
  int num_rotations = 0;

  while (max_elem > eps)
  {
    Jacobi_rot(A, max_elem);
    num_rotations += 1;
    //printf("max_elem outside function: %g\n", max_elem);
  }

  vec eigvals = sort(A.diag());

  ofstream myfile;
  char *project2_e_data;
  myfile.open ("project2_e_data.txt");
  for (int i=0; i < 4; i++)
  {
    myfile << eigvals(i) << endl;
  }

  myfile.close();
  printf("Number of rotations performed: %i\n", num_rotations );
  return 0;
}
