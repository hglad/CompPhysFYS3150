#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

class Planet
{
  private:
    double pi = M_PI;
    double GM = 4*pi*pi;
    double M_Sun = 1;

  public:
    mat data_mats(int n, double x_0, double y_0, double vx_0, double vy_0, double M)
    {
      mat s = zeros(n, 2);
      mat v = zeros(n, 2);
      mat a = zeros(n, 2);

      s(0,0) = x_0;
      s(0,1) = y_0;
      v(0,0) = vx_0;
      v(0,1) = vy_0;

      double r_0 = sqrt(s(0,0)*s(0,0) + s(0,1)*s(0,1));
      double F_0 = GM*M/(r_0*r_0);
      a(0,0) = -F_0/M * s(0,0);
      a(0,1) = -F_0/M * s(0,1);


    }

};
