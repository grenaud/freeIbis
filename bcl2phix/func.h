#include <string>
#include <iostream>


using namespace std;

inline char complement(const char c){
    if(c ==    'A')
	return 'T';

    if(c ==    'C')
	return 'G';

    if(c ==    'G')
	return 'C';

    if(c ==    'T')
	return 'A';



    if(c ==    'a')
	return 't';

    if(c ==    'c')
	return 'g';

    if(c ==    'g')
	return 'c';

    if(c ==    't')
	return 'a';



    if(c ==    'N')
	return 'N';

    cerr<<"bcl2phix: complement: Invalid base pair="<<c<<endl;
    exit(1);
}


inline string reverseComplement(const string & inputString){
    string toReturn="";
    if(inputString.size() >0 )	
	for(int i=(inputString.size()-1);i>=0;i--){
	    toReturn+=complement( inputString[i] );
	}

    return toReturn;
}
