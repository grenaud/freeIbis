#ifndef INCLUDED_INPUT_H
#define INCLUDED_INPUT_H

#include <math.h>
#include <string>

extern int verbosity;

// Generalized input, modelled to mimic the Firecrest-style files.  
// We run an outer loop over clusters (characterized by coordinates)
// and an inner loop over cycles.  
class Input
{
    public:
        virtual ~Input() {}

        //! \brief moves to next cluster and reads coordinates
        //!
        //! If another cluster exists, reads its coordinates into
        //! arguments and returns true.  Else returns false and
        //! arguments become undefined.  Behaviour is undefined if a
        //! previous next_cluster() did not return false.
        virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) = 0 ;

        //! \brief moves to next cycle and reads four intensities.
        //!
        //! If another cycle exists, this reads the four intensity
        //! values and returns true.  Else return false and the in/out
        //! parameters become undefined.  The first call after
        //! next_cluster() reads the first cycle.   
        virtual bool next_cycle( double* pa, double* pc, double* pg, double* pt ) = 0 ;
} ;

// Generalized input of coordinates only, for formats that have separate
// position files (which exist in three variants now).
class PosnInput 
{
    public:
        virtual ~PosnInput() ;
        virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) = 0 ;
        virtual int lane() const = 0 ;
        virtual int tile() const = 0 ;
} ;

PosnInput *new_ipar_posn_input( const std::string& index_fn, int coordtype ) ;
PosnInput *new_locs_input( const std::string& index_fn, int coordtype ) ;
PosnInput *new_clocs_input( const std::string& index_fn ) ;

// Construct correct position reader depending on file type.  As much as
// it sucks, the only way to identify the file type is the name. :(
PosnInput* new_posn_input( const std::string& index_fn, int coordtype ) ;

Input *new_firecrest_input( size_t nmodels, const char* fp ) ;
Input *new_ipar_input( size_t nmodels, const char* index, int coordtype, const char* fp ) ;
Input *new_cif_input( size_t nmodels, PosnInput* pinput, const char* fp ) ;

void split_filename( const std::string& fn, int nfields, int *plane, int *ptile ) ;

static inline int do_round( int coordtype_, float x )
{
    switch( coordtype_ )
    {
        case 1: return int( round( x ) ) ;
        case 2: return int( floor( fabs( x ) ) ) ;
        case 3: return int( round( 10 * x + 1000 ) ) ;
        default: throw "unknown coordinate format" ;
    }
}

#endif

