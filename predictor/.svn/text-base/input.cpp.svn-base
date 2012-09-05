#include "input.h"

#include <iostream>
#include <sstream>
#include <string>

int verbosity = 0 ;

PosnInput::~PosnInput() {}

static bool ends_with( const std::string& s, const std::string& t )
{
    return s.substr( s.length() - t.length() ) == t ;
}

PosnInput* new_posn_input( const std::string& fn, int ctype )
{
    if( ends_with( fn, "_pos.txt" ) ) return new_ipar_posn_input( fn, ctype ) ;
    if( ends_with( fn, ".locs" ) ) return new_locs_input( fn, ctype ) ;
    if( ends_with( fn, ".clocs" ) ) return new_clocs_input( fn ) ;
    throw "cannot infer format of position file" ;
}

void split_filename( const std::string& fn, int nfields, int *plane, int *ptile )
{

    // New index file does not contain lane and tile, try to extract
    // from filename.  [Okay, granted, the split function was ugly, but
    // so are the C++ string algorithms...]
    std::string::size_type p = std::string::npos+1 ;
    for( int i = 0 ; i != nfields && std::string::npos != (p = fn.rfind( '_', p-1 )) ; ++i ) ;

    std::istringstream coords( fn.substr( p ) ) ;
    char junk ;
    if( coords >> junk >> *plane >> junk >> *ptile ) 	{
        if(verbosity>=1) 
            std::cout << "Reconstructed lane and tile information from "
                "filename: " << *plane << ',' << *ptile << '.' << std::endl ;
    } else{
        std::cout << "Reconstruction of lane and tile from '" << fn.substr(p) << "' failed, will use 0,0." << std::endl ;
        *plane = 0 ;
        *ptile = 0 ;
    }
}

