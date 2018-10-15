#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

class planet
{
  private:
    double pi = M_PI;
    double GM = 4*pi*pi;
    double M_Sun = 1;

  public:
    // Generate matrix for planet number j
    void data_mats(int j, double x_0, double y_0, double vx_0, double vy_0, double M)
    {
      // Store initial values in matrices
      s(0,0,j) = x_0;
      s(0,1,j) = y_0;
      v(0,0,j) = vx_0;
      v(0,1,j) = vy_0;

      double r_0 = sqrt(s(0,0,j)*s(0,0,j) + s(0,1,j)*s(0,1,j));
      double F_0 = GM*M/(r_0*r_0);
      a(0,0,j) = -F_0/M * s(0,0,j)/r_0;
      a(0,1,j) = -F_0/M * s(0,1,j)/r_0;

      return;

    }

};
