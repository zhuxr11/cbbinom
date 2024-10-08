#ifndef CBBINOM_ZEROIN
#define CBBINOM_ZEROIN

#include <float.h>
#include <math.h>

#include <R_ext/Applic.h>

namespace cbbinom
{

#ifndef EPSILON
#define EPSILON DBL_EPSILON
#endif

class UnirootEqn
{
public:
  virtual double operator () (const double& x) const = 0;
};

inline double cpp_uniroot(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    const UnirootEqn* const f,	/* Function under investigation	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit				/* Max # of iterations */
) {
  double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
  double tol;
  int maxit;

  a = ax;  b = bx;
  c = a;   fc = fa;
  maxit = *Maxit + 1; tol = * Tol;

  /* First test if we have found a root at an endpoint */
  if(fa == 0.0) {
    *Tol = 0.0;
    *Maxit = 0;
    return a;
  }
  if(fb ==  0.0) {
    *Tol = 0.0;
    *Maxit = 0;
    return b;
  }

  while(maxit--)		/* Main iteration loop	*/
  {
    double prev_step = b-a;		/* Distance from the last but one
   to the last approximation	*/
  double tol_act;			/* Actual tolerance		*/
  double p;			/* Interpolation step is calcu- */
  double q;			/* lated in the form p/q; divi-
   * sion operations is delayed
   * until the last moment	*/
  double new_step;		/* Step at this iteration	*/

  if( fabs(fc) < fabs(fb) )
  {				/* Swap data for b to be the	*/
  a = b;  b = c;  c = a;	/* best approximation		*/
  fa=fb;  fb=fc;  fc=fa;
  }
  tol_act = 2*EPSILON*fabs(b) + tol/2;
  new_step = (c-b)/2;

  if( fabs(new_step) <= tol_act || fb == (double)0 )
  {
    *Maxit -= maxit;
    *Tol = fabs(c-b);
    return b;			/* Acceptable approx. is found	*/
  }

  /* Decide if the interpolation can be tried	*/
  if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
  && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
   * Interpolation may be tried	*/
  double t1,cb,t2;
    cb = c-b;
    if( a==c ) {		/* If we have only two distinct	*/
  /* points linear interpolation	*/
  t1 = fb/fa;		/* can only be applied		*/
  p = cb*t1;
  q = 1.0 - t1;
    }
    else {			/* Quadric inverse interpolation*/

  q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
  q = (q-1.0) * (t1-1.0) * (t2-1.0);
    }
    if( p>(double)0 )		/* p was calculated with the */
  q = -q;			/* opposite sign; make p positive */
  else			/* and assign possible minus to	*/
  p = -p;			/* q				*/

  if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
  && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
  new_step = p/q;			/* it is accepted
   * If p/q is too large then the
   * bisection procedure can
   * reduce [b,c] range to more
   * extent */
  }

  if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
  if( new_step > (double)0 )	/* than tolerance		*/
  new_step = tol_act;
  else
    new_step = -tol_act;
  }
  a = b;	fa = fb;			/* Save the previous approx. */
  b += new_step;	fb = (*f)(b);	/* Do step to a new approxim. */
  if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
    /* Adjust c for it to have a sign opposite to that of b */
    c = a;  fc = fa;
  }

  }
  /* failed! */
  *Tol = fabs(c-b);
  *Maxit = -1;
  return b;
}

} // namespace cbbinom

#endif // CBBINOM_ZEROIN
