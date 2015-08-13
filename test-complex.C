

/* check that exp(i*pi) == -1 */
#include <math.h>   /* for atan */
#include <complex>
using namespace std;
main() 
{
  double pi = 4*std::atan(1.0);
  complex<double> I(0,1);
  complex<double> z = exp(I*pi);
  printf("%f+%f*i\n", z.real(), z.imag());//creal(z), cimag(z));
}
