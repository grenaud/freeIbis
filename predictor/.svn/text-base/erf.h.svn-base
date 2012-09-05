#ifndef INCLUDED_ERF_H
#define INCLUDED_ERF_H

#include <math.h>

static inline void add_logs( double &x, double y ) 
{
    if( x > y ) x += log1p( exp( y-x ) ) ;
    else x = y + log1p( exp( x-y ) ) ;
}

double log_cum_density(double x,double mu, double sigma) ;

#endif
