#include "input.h"

#include "../gzstream/gzstream.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std ;

class IparPosnInput : public PosnInput
{
    private:
    int lane_, tile_ ;
    ifstream index_ ;
    int coordtype_ ;

public:
    IparPosnInput( const string& index_fn, int coordtype ) ;
    virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) ;
    virtual int lane() const { return lane_ ; }
    virtual int tile() const { return tile_ ; }
} ;

class IparInput : public Input
{
private:
    IparPosnInput pinput_ ;
    size_t ncycles_, nclusters_ ;
    size_t cur_cluster_, cur_cycle_ ;

    // compact representation of intensities, same order as we read
    // them in, that is 
    // for each cyle
    //   for each cluster
    //     for each channel
    //       one float
    vector< float > intensities_ ;

public:
    IparInput( size_t nmodels, const char* index, int coordtype, const char* fp ) ;
    virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) ;
    virtual bool next_cycle( double* pa, double* pc, double* pg, double* pt ) ;
} ;


PosnInput *new_ipar_posn_input( const string& index_fn, int coordtype )
{
    return new IparPosnInput( index_fn, coordtype ) ;
}

IparPosnInput::IparPosnInput( const string& index_fn, int coordtype ) :
	index_( index_fn.c_str() ), coordtype_( coordtype )
{
	if(verbosity>=1) cout << "Opening IPAR/CIF index file: " << index_fn << endl;
	split_filename( index_fn,3, &lane_, &tile_ ) ;
}

bool IparPosnInput::next_cluster( int *plane, int *ptile, int *px, int* py )
{
	string line ;
	if( !getline( index_, line ) ) return false ;

	// first try: old index format with four entries (_idx.txt file)
	istringstream ss( line ) ;
	if( ss >> *plane >> *ptile >> *px >> *py ) return true ;

	// failed: might be new index with two float coords (_pos.txt)
	ss.str( line ) ;
	ss.clear() ;
	float x, y ;
	if( ss >> x >> y ) 
	{
		*plane = lane_ ;
		*ptile = tile_ ;
        *px = do_round( coordtype_, x ) ;
        *py = do_round( coordtype_, y ) ;
		return true ;
	}

	// fall through: must be end of file (or an error condition I don't
	// want to think about)
	return false ;
}

Input *new_ipar_input( size_t nmodels, const char* index, int coordtype, const char* fp )
{
    return new IparInput( nmodels, index, coordtype, fp ) ;
}

IparInput::IparInput( size_t nmodels, const char* index, int coordtype, const char* fp ) :
    pinput_( index, coordtype ), ncycles_(0), nclusters_(0), cur_cluster_(0)
{
    if(verbosity>=1) { cout << "Reading IPAR file to memory: " << fp << endl; }

    igzstream ints( fp ) ;
    size_t cluster = 0 ;

    // one header ('#')
    //  one line per cluster
    // one footer line per cycle ('#')
    string line ;
    if( !getline( ints, line ) || line.empty() || line[0] != '#' )
	throw "intensities files has no header" ;

    while( getline( ints, line ) && ncycles_ != nmodels )
	{
	    if( !line.empty() && line[0] == '#' )
		{
		    if( nclusters_ )
			{
			    if( cluster != nclusters_ )
				throw "inconsistent number of clusters in IPAR file" ;
			}
		    else nclusters_ = cluster ;
		    cluster = 0 ;
		    ++ncycles_ ;
		}
	    else
		{
		    istringstream ss( line ) ;
		    float a, c, g, t ;
		    if( ss >> a >> c >> g >> t ) 
			{
			    intensities_.push_back( a ) ;
			    intensities_.push_back( c ) ;
			    intensities_.push_back( g ) ;
			    intensities_.push_back( t ) ;
			    ++cluster ;
			}
		    else throw "parse error in IPAR file" ;
		}
	}

    if( verbosity ) 
	cout << "Got " << nclusters_ << " clusters of " << ncycles_ << " cycles each in "
	     << (intensities_.size() * sizeof(float) >> 20) << "MiB" << endl ;
}

bool IparInput::next_cluster( int *plane, int *ptile, int *px, int* py )
{
    ++cur_cluster_ ;
    cur_cycle_ = 0 ;
    return pinput_.next_cluster( plane, ptile, px, py ) ;
}

bool IparInput::next_cycle( double* pa, double* pc, double* pg, double* pt )
{
    if( cur_cycle_ == ncycles_ ) return false ;
    *pa = intensities_[ 0 + 4 * ( cur_cluster_-1 + nclusters_ * cur_cycle_ ) ] ;
    *pc = intensities_[ 1 + 4 * ( cur_cluster_-1 + nclusters_ * cur_cycle_ ) ] ;
    *pg = intensities_[ 2 + 4 * ( cur_cluster_-1 + nclusters_ * cur_cycle_ ) ] ;
    *pt = intensities_[ 3 + 4 * ( cur_cluster_-1 + nclusters_ * cur_cycle_ ) ] ;
    ++cur_cycle_ ;
    return true ;
}
