#include <math.h>
#include <stdio.h>


/** Square root of Pi = 3.14159
* Would for example be defined in gsl/gsl_math.h .
*/
#define M_SQRTPI 1.77245385090551602729816748334
/*************************************
* @brief Compute Boys’ function F_m(z) via the asymptotic high-z expansion.
* @param m The non-negative integer index.
* The function supports only the range 0<=m<=128.
* @param z The non-negative radial distance
* @param N The number of terms in the 1/z expansion.
* @return F_m(z) with approximately double precision accuracy.
* @author Richard J. Mathar
* @since 2017-09-17 Richard J. Mathar
* https://vixra.org/pdf/1709.0304v1.pdf
*/
static double boysHighz(int m, double z, int N) {
	double a = m + 0.5;
	double f = 1.0;
	double fn = 1.0;
	double gam = M_SQRTPI;
	int n = 1;
	for (; n <= N; n++) {
		fn *= (a - n) / z;
		f += fn;
	}
	f *= -a * exp(-z) / z;
	/* add Gamma(1+a)/z^a, where Gamma(1/2)=sqrt(Pi)
	*/
	int j = 0;
	for (; j <= m; j++)
		gam *= j + 0.5;
	f += gam / pow(z, a);
	/* divide through 2m+1 to convert 1F1(a,a+1,-z) to F_m(z) */
	return f / (2 * a);
}
#undef M_SQRTPI
/**********************************
* @brief Compute Boys’ function F_m(z) via Gauss-Jacobi quadrature.
* @param m The non-negative integer index.
* The function supports only the range 0<=m<=128.
* @param z The non-negative radial distance
* @return F_m(z) with approximately double precision accuracy.
* If m or z are outside the range of supported values, 0 is returned.
* @author Richard J. Mathar
* @since 2017-09-17 Richard J. Mathar
*/
double boysGaussJacobi(int m, double z) {
#include "boysGaussJacobiData.h"
	/* if gj is not assigned here, either m is not
	* in the supported range or z is outside the region of convergence
	*/
	if (gj && N > 0) {
		double F = 0.0;
		/* N abscissae are in gj[0,2,4...] and N weights in gj[1,3,5....]
		*/
		int i = 0;
		for (; i < N; i++)
			F += gj[1 + 2 * i] * exp(-z * gj[2 * i]);
		/* gather the factor 1/2 left over from the u->v substitution
		*/
		return F / 2.0;
	}
	else if (N > 0) {
		/* call the high-z expansion with N terms of the Laurent series */
		return boysHighz(m, z, N);
	}
	else {
		/* large m: Kummer transformation assuming z/m is small.
		* 1F1(a,a+1,-z) = exp(-z)* 1F1(1;a+1;z) */
		double a = m + 0.5;
		double f = 1.0;
		double fn = 1.0;
		int n = 1;
		for (; ; n++)
		{
			fn *= z / (a + n);
			f += fn;
			if (fabs(fn / f) < 1.1e-16)
				break;
		}
		/* Kummer’s factor exp(z). Factor 1/(2m+1) to get F_m
		*/
		return f * exp(-z) / (2 * a);
	}
} /* boysGaussJacobi */