#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

class planet
{
  private:
    double pi = M_PI;
    double GM = 4*pi*pi;

  public:
    // Needs vectors to hold initial values since we consider multiple planets
    void data_mats(int n, mat& s, mat& v, mat& a, vec s_0, vec v_0, double M)
    {
      // Loop through all planets
      s(0, 0) = s_0(0);
      s(0, 1) = s_0(1);
      v(0, 0) = v_0(0);
      v(0, 1) = v_0(1);

      double r_0 = sqrt(s(0,0)*s(0,0) + s(0,1)*s(0,1));
      double F_0 = GM*M/(r_0*r_0);
      a(0,0) = -F_0/M * s(0,0)/r_0;
      a(0,1) = -F_0/M * s(0,1)/r_0;

      return;

    }

};
