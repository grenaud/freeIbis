#include "erf.h"

// Log of cumulative density of normal distribution.
// Originally, we calculated the cumulative density as
//
//   1/2 + 1/2 * erf( (x-mu)/(sigma*sqrt2pi) ) 
//
// but that runs into trouble for very small x.  Since we take a
// logarithm soon afterwards anyway, it's better to directly calculate
// the logarithm in a way that doesn't cause loss.  (Note: natural log,
// not Phred scale!)

// Code liberated from GSL (licensed under GPL).
// See Hart et al, Computer Approximations, John Wiley and Sons, New York (1968)

static const double sqrt2pi = sqrt( 2*M_PI ) ;

static double log_erfc8_sum(double x)
{
    // estimates erfc(x) valid for 8 < x < 100
    // This is based on index 5725 in Hart et al

    static double P[] = {
	2.97886562639399288862,
	7.409740605964741794425,
	6.1602098531096305440906,
	5.019049726784267463450058,
	1.275366644729965952479585264,
	0.5641895835477550741253201704
    };
    static double Q[] = {
	3.3690752069827527677,
	9.608965327192787870698,
	17.08144074746600431571095,
	12.0489519278551290360340491,
	9.396034016235054150430579648,
	2.260528520767326969591866945,
	1.0
    };
    double num=0.0, den=0.0;
    int i;

    num = log( P[5] ) ;
    for (i=4; i>=0; --i) {
	num += log(x) ;
	add_logs( num, log( P[i] ) ) ;
    }
    den = log( Q[6] ) ;
    for (i=5; i>=0; --i) {
	den += log(x) ;
	add_logs( den, log( Q[i] ) ) ;
    }

    return num - den;
}

static double lerfc( double x )
{
    return x > 8.0
	? log_erfc8_sum( x ) - x*x 
	: log( erfc( x ) ) ;
}

double log_cum_density(double x,double mu, double sigma)
{
    return lerfc( (mu-x)/(sigma*sqrt2pi) ) - log(2) ;
}

