#include "input.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include <dirent.h>
#include <fcntl.h>
#include <fnmatch.h>
#include <stdint.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std ;

static inline int decode_int( int size, const uint8_t *base )
{
    int r = 0, s = 0 ;
    while( size > 1 )
	{
	    r |= int( uint8_t( *base ) ) << s ;
	    s += 8 ;
	    --size ;
	    ++base ;
	}
    r |= int( int8_t( *base ) ) << s ;
    return r ;
}

// Input from CIF files.  These days, CIF files can be combined with
// textual position files (GA IIx), clogs files (HiSeq) or locs files
// (MiSeq).  We delegate to an instance of PosnInput for the details.
class CifInput : public Input
{
private:
	        PosnInput *pinput_ ;

		// set of CIF files, mapped into memory
		std::vector< const uint8_t* > blob_ ;
		std::vector< int > size_ ;
		unsigned nclusters_, cur_cluster_, cur_cycle_, data_size_ ;

	public:
		CifInput( size_t nmodels, PosnInput*, const char* fp ) ;
   	        ~CifInput() ;
		virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) ;
		virtual bool next_cycle( double* pa, double* pc, double* pg, double* pt ) ;
} ;

Input *new_cif_input( size_t nmodels, PosnInput* pinput, const char* fp )
{
    return new CifInput( nmodels, pinput, fp ) ;
}

CifInput::CifInput( size_t nmodels, PosnInput* pinput, const char* fp ) :
    pinput_( pinput ), nclusters_(0), cur_cluster_(0), cur_cycle_(0)
{
    if(verbosity>=1) {
	cout << "Reading CIF files to memory: " << fp << "/C*.1/s_"
	     << pinput_->lane() << '_' << pinput_->tile() << ".cif" << endl ;
    }
	
    long mem = 0 ;
	
    // what to do?
    // - open directory, scan it
    // - descend into directory, read correct file
    // - first read, later verify header information
    DIR *dir = opendir( fp ) ;
    if( !dir ) throw "could not open CIF directory" ;
    while( struct dirent *de = readdir( dir ) )
	{
	    if( 0 == fnmatch( "C*.1", de->d_name, 0 ) )
		{
		    int cycleno ;
		    sscanf( de->d_name, "C%d.1", &cycleno ) ;
		    if( (unsigned)cycleno > nmodels ) continue ; // limit to tnumber of available models
	
		    if(verbosity>=2) cout << cycleno << ' ' << flush ;
		    --cycleno ; // stupid idiots start counting at 1
	
		    if( blob_.size() <= (unsigned)cycleno ) { blob_.resize( cycleno+1 ) ; size_.resize( cycleno+1 ) ; }
		    ostringstream fname ;
		    fname << fp << '/' << de->d_name << "/s_" << pinput_->lane() << '_' << pinput_->tile() << ".cif" ;
		    int fd = open( fname.str().c_str(), O_RDONLY ) ;
		    if( fd == -1 && errno != ENOENT ) throw fname.str() + ": " + strerror( errno ) ;
	
		    struct stat the_stat ;
	            void *p = 0;
	            try
			{
			    if( fd == -1 ) throw fname.str() + ": " + strerror( errno ) ;
			    fstat( fd, &the_stat ) ;
			    if( the_stat.st_size < 13 ) throw fname.str() + ": file too small" ;
	
			    void *p = mmap( 0, the_stat.st_size, PROT_READ, MAP_SHARED, fd, 0 ) ;
			    if( p == (void*)-1 ) throw "error when memory-mapping " + fname.str() ;
			    close( fd ) ;
	
			    uint8_t* raw = (uint8_t*)p ;
			    blob_[ cycleno ] = (const uint8_t*)raw ;
			    size_[ cycleno ] = the_stat.st_size ;
			    mem += the_stat.st_size ;
	
			    // Check header information: what do we do if it is wrong?
			    if( raw[0] != 'C' || raw[1] != 'I' || raw[2] != 'F' || raw[3] != 1 )
				throw fname.str() + ": not a CIF/1 file" ;
	
			    if( raw[4] != 1 && raw[4] != 2 && raw[4] != 4 )
				throw fname.str() + ": this data size is not allowed" ;
	
			    if( decode_int( 2, raw+5 ) != cycleno+1 || decode_int( 2,raw+7 ) != 1 )
				cerr << "warning: " << fname.str()
				     << ": does not appear to contain the correct cycle, continuing anyway." ;
	
			    if( nclusters_ )
				{
				    if( raw[4] != data_size_ || decode_int( 4, raw+9 ) != (int)nclusters_ )
					throw fname.str() + ": header information is inconsistent with sibling files" ;
				}
			    else
				{
				    data_size_ = raw[4] ;
				    nclusters_ = decode_int( 4, raw+9 ) ;
				}
	
	                    if( the_stat.st_size < 13 + 4 * data_size_ * nclusters_ )
				throw  fname.str() + ": file too small" ;
	
			    // Advise memory manager that we're going to read four
			    // separate streams from the file.  This is supposed to
			    // improve performance somewhat.
			    int psize = sysconf( _SC_PAGE_SIZE ) ;
			    int stride = (data_size_ * nclusters_ + psize - 1) / psize * psize ;
			    madvise( raw,            stride, MADV_SEQUENTIAL ) ;
			    madvise( raw +   stride, stride, MADV_SEQUENTIAL ) ;
			    madvise( raw + 2*stride, stride, MADV_SEQUENTIAL ) ;
			    madvise( raw + 3*stride, stride, MADV_SEQUENTIAL ) ;
			}
	            catch( const std::string& e )
			{
			    std::cerr << e << ", will substitute zero intensities" << std::endl ;
			    if( p )
				{
				    munmap( p, the_stat.st_size ) ;
				    mem += the_stat.st_size ;
				}
			    blob_[ cycleno ] = 0 ;
			}
		}
	}
    if(verbosity>=2) cout << endl ;
    closedir( dir ) ;
	
    for( unsigned i = 0 ; i != blob_.size() ; ++i )
	{
	    if( !blob_[i] )
		std::cerr << "couldn't find CIF file for cycle " << i << ", will assume zero intesities." << std::endl ;
	}
	
    if( verbosity )
	cout << "Got " << nclusters_ << " clusters with " << blob_.size()
	     << " cycles, data size " << data_size_ << ", " << (mem >> 20) << "MiB" << endl;
}

CifInput::~CifInput()
{
	for( unsigned i = 0 ; i != blob_.size() ; ++i )
        if( blob_[i] )
            munmap( (void*)blob_[i], size_[i] ) ;
        delete pinput_ ;
}

bool CifInput::next_cluster( int *plane, int *ptile, int *px, int* py )
{
    ++cur_cluster_ ;
    cur_cycle_ = 0 ;
    return pinput_->next_cluster( plane, ptile, px, py ) ;
}

bool CifInput::next_cycle( double* pa, double* pc, double* pg, double* pt )
{
    if( cur_cycle_ == blob_.size() ) return false ;

    int off = data_size_ * (cur_cluster_-1) + 13 ;
    int stride = data_size_ * nclusters_ ;

    if( blob_[ cur_cycle_ ] )
    {
        *pa = decode_int( data_size_, &blob_[ cur_cycle_ ][ off              ] ) ;
        *pc = decode_int( data_size_, &blob_[ cur_cycle_ ][ off +     stride ] ) ;
        *pg = decode_int( data_size_, &blob_[ cur_cycle_ ][ off + 2 * stride ] ) ;
        *pt = decode_int( data_size_, &blob_[ cur_cycle_ ][ off + 3 * stride ] ) ;

        int m = std::numeric_limits<short>::min() ;
        if( *pa == m || *pc == m || *pg == m || *pt == m )
            *pa = *pc = *pg = *pt = 0 ;
    }
    else *pa = *pc = *pg = *pt = 0 ;

    ++cur_cycle_ ;
    return true ;
}

