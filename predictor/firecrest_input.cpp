#include "input.h"

#include "../gzstream/gzstream.h"

#include <iostream>
#include <sstream>

using namespace std ;

// Firecrest-style input.  This is a single file and the iterator
// interface reduces to incremental reading of the input.
class FirecrestInput : public Input
{
    private:
        igzstream intensities_ ;
        istringstream cur_line_ ;
        size_t nmodels_, cycle_ ;

    public:
        FirecrestInput( size_t nmodels, const char* fp ) ;
        virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) ;
        virtual bool next_cycle( double* pa, double* pc, double* pg, double* pt ) ;
} ;


FirecrestInput::FirecrestInput( size_t nmodels, const char* fp ) : intensities_( fp ), nmodels_(nmodels), cycle_(0)
{
    if(verbosity>=1) cout << "Reading Firecrest file: " << fp << endl ;
}

bool FirecrestInput::next_cluster( int *plane, int *ptile, int *px, int* py ) 
{
    string line ;
    if( !getline( intensities_, line ) ) return false ;
    cur_line_.str( line ) ;
    cur_line_.clear() ;
    cycle_ = 0 ;
    return cur_line_ >> *plane >> *ptile >> *px >> *py ; 
}

bool FirecrestInput::next_cycle( double* pa, double* pc, double* pg, double* pt )
{
    if( cycle_ == nmodels_ ) return false ;
    ++cycle_ ;
    return cur_line_ >> *pa >> *pc >> *pg >> *pt ;
}

Input *new_firecrest_input( size_t nmodels, const char* fp )
{
    return new FirecrestInput( nmodels, fp ) ;
}

