#include <math.h>

float gasdev(long *idum)
// Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum) as the source of uniform deviates.

{
    float ran1(long *idum);
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if (*idum < 0) iset=0;          // Reinitialize.	      
    if (iset == 0) {                // We don’t have an extra deviate handy, so
      do {
	v1=2.0*ran1(idum)-1.0;  // pick two uniform numbers in the square extending from -1 to +1 in each direction,
	v2=2.0*ran1(idum)-1.0; // see if they are in the unit circle, 
	rsq=v1*v1+v2*v2;     // and if they are not, try again
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      // Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
      
      gset=v1*fac;
      iset=1;      // Set flag.
      return (float)(v2*fac);
    } else {         // We have an extra deviate handy,
      iset=0;   // so unset the flag,
      return (float)gset;   // and return it.
    }
}
