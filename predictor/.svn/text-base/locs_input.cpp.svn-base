#include "input.h"

#include <fstream>
#include <iostream>
#include <string>
#include <stdint.h>

using namespace std ;

/* Input from "locs" file.  Python code for comparison:

def read_locs(filename):
  global coordconv
  infile = open(filename,'rb')
  infile.read(8) # First 8 Byte are unused
  clusters = to_int(infile.read(4))
  binvalues = array.array('f')
  binvalues.read(infile, clusters * 2)
  data = N.array(binvalues, typecode=N.Float)
  data = N.reshape(data, ( clusters ,2 ))
  for x,y in data:
    yield coordconv(x),coordconv(y)
  infile.close()
  raise StopIteration

*/

// Note that due to lack of documentation, this code works only if on
// the current platform, a float is a 32bit IEEE floating point number.
class LocsInput : public PosnInput
{
    private:
		int lane_, tile_ ;
        ifstream file_ ;
        int nclusters_ ;
        int coordtype_ ;

	public:
		LocsInput( const string& index_fn, int coordtype ) ;
		virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) ;
        virtual int lane() const { return lane_ ; }
        virtual int tile() const { return tile_ ; }
} ;

PosnInput *new_locs_input( const string& index_fn, int coordtype )
{
    return new LocsInput( index_fn, coordtype ) ;
}

LocsInput::LocsInput( const string& index_fn, int coordtype )
    : file_( index_fn.c_str() ), coordtype_( coordtype )
{
	if(verbosity>=1) cout << "Opening LOCS index file: " << index_fn << endl;
	split_filename( index_fn, 2,&lane_, &tile_ ) ;

    uint32_t gunk[3] ;
    file_.read( (char*)gunk, 12 ) ;
    nclusters_ = gunk[2] ;
}

bool LocsInput::next_cluster( int *plane, int *ptile, int *px, int* py )
{
    if( !nclusters_ ) return false ;
    *plane = lane_ ;
    *ptile = tile_ ;

    float cc[2] ;
    file_.read( (char*)cc, 8 ) ;
    *px = do_round( coordtype_, cc[0] ) ;
    *py = do_round( coordtype_, cc[1] ) ;
    return file_ ;
}

/* Input from "clocs" compressed locs file.  Python code for comparison:

def read_clocs(filename):
  EXPECTED_CLOCS_VERSION = 1
  BLOCK_SIZE = 25
  IMAGE_WIDTH = 2048
  BLOCKS_PER_LINE = (IMAGE_WIDTH + BLOCK_SIZE - 1) / BLOCK_SIZE
  totalBlocks = 0
  currentBlock = 0
  currentBlockUnreadClusters = 0

  infile = open(filename,'rb')
  clocsVersion = ord(infile.read(1))
  totalBlocks = to_int(infile.read(4))
  currentBlockUnreadClusters = ord(infile.read(1))
  currentBlock+=1

  while (currentBlock < totalBlocks or ( currentBlock == totalBlocks and currentBlockUnreadClusters > 0)):
     while (currentBlockUnreadClusters == 0 and currentBlock < totalBlocks):
        currentBlockUnreadClusters = ord(infile.read(1))
        currentBlock += 1
     dx = ord(infile.read(1))
     dy = ord(infile.read(1))
     x = 10 * BLOCK_SIZE * ((currentBlock - 1) % BLOCKS_PER_LINE) + dx + 1000;
     y = 10 * BLOCK_SIZE * ((currentBlock - 1) / BLOCKS_PER_LINE) + dy + 1000;
     yield x,y
     currentBlockUnreadClusters -= 1
  infile.close()
  raise StopIteration

*/

class ClocsInput : public PosnInput
{
    private:
		int lane_, tile_ ;
        ifstream file_ ;
        int nblocks_ ;
        int blocknum_ ;
        unsigned char nclusters_ ;

	public:
		ClocsInput( const string& index_fn ) ;
		virtual bool next_cluster( int *plane, int *ptile, int *px, int* py ) ;
        virtual int lane() const { return lane_ ; }
        virtual int tile() const { return tile_ ; }
} ;

PosnInput *new_clocs_input( const string& index_fn )
{
    return new ClocsInput( index_fn ) ;
}

ClocsInput::ClocsInput( const string& index_fn )
    : file_( index_fn.c_str() )
{
	if(verbosity>=1) cout << "Opening CLOCS index file: " << index_fn << endl;
	split_filename( index_fn, 2,&lane_, &tile_ ) ;

    char gunk[5] ;
    file_.read( gunk, 5 ) ;
    switch( gunk[0] )
    {
        case 1: 
            nblocks_ = (int)gunk[1] | (int)gunk[2] << 8 | (int)gunk[3] << 16 | (int)gunk[4] << 24 ;
            nclusters_ = 0 ;
            blocknum_ = -1 ;
            break ;
        default:
            throw "CLOCS file has unknown version number" ;
    }
}

bool ClocsInput::next_cluster( int *plane, int *ptile, int *px, int* py )
{
    static const int block_size = 25 ;
    static const int image_width = 2048 ;
    static const int blocks_per_line = (image_width + block_size - 1) / block_size ;

    while( !nclusters_ ) {
        if( !nblocks_ ) return false ;
        file_.read( (char*)&nclusters_, 1 ) ;
        --nblocks_ ;
        ++blocknum_ ;
    }
    nclusters_--;

    unsigned char dd[2] ;
    file_.read( (char*)&dd, 2 ) ;

    *plane = lane_ ;
    *ptile = tile_ ;
    *px = 10 * block_size * (blocknum_ % blocks_per_line) + dd[0] + 1000 ;
    *py = 10 * block_size * (blocknum_ / blocks_per_line) + dd[1] + 1000 ;
    return file_ ;
}


