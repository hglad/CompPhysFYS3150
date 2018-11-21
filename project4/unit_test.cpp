#define CATCH_CONFIG_MAIN
#include "ising.h"
#include "catch.hpp"

/* Test that the calculated expectation values of a 2 x 2 lattice match the
analytical values
compile with:
g++ unit_test.cpp ising.cpp -o test.x -DARMA_DONT_USE_WRAPPER -lblas -llapack
*/

TEST_CASE("Check if numerical solution for a 2x2 at T = 1 match analytical values")
{
  // Analytical
  vec analytical = zeros(5);

  analytical[0] = -7.98393;
  analytical[1] = 3.99464;
  analytical[2] = 15.97322;
  analytical[3] = 0.12833;
  analytical[4] = 15.97322;

  vec ValueSums = zeros(5);
  int numMC = 500000;
  int L = 2; int T = 1;

  mt19937_64 gen(5);
  mat S = init_spins(L, gen, 0);
  map<double, double> w = transitions(T);

  double energy, magmom;
  int counter = 0;
  init_params(S, energy, magmom);   // initial energy, magnetic momentum
  for (int k = 1; k < numMC; k++)
  {
    MC_cycle(S, L, counter, energy, magmom, w, gen);

    ValueSums(0) += energy; ValueSums(1) += energy*energy;
    ValueSums(2) += magmom; ValueSums(3) += magmom*magmom;
    ValueSums(4) += fabs(magmom);
  }
  double E = ValueSums(0)/numMC;
  double absM = ValueSums(4)/numMC;
  double M = ValueSums(2)/numMC;
  double M2 = ValueSums(3)/numMC;
  double C_V = (ValueSums(1)/numMC - E*E)/(T*T);
  double chi = (M2 - M*M)/T;

  vec means = {E, absM, M2, C_V, chi};

  for (int i=0; i < 5; i++)
  {
    REQUIRE( means(i) == Approx(analytical(i)).epsilon(0.03) );
  }
}
