/*
 * BAMWriter
 * Date: Sep-04-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef BAMWriter_h
#define BAMWriter_h

#include <limits.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <inttypes.h>

using namespace std;

//CONSTANTS
const unsigned char zero = 0;
const uint16_t flagSingleReads =  4; // 00000100
const uint16_t flagFirstPair   = 77; // 01001101
const uint16_t flagSecondPair  =141; // 10001101

const unsigned char valueN=15;
const unsigned char valueA=1;
const unsigned char valueC=2;
const unsigned char valueG=4;
const unsigned char valueT=8;


class BAMWriter{
private:
    ofstream * myfile;
    unsigned char bp2number(char b);
    void writeZero();
    void writeTag(string & tag);
    inline void writeSequence(uint16_t flagToUse,string & name,string & firstSeq,string & firstQual,string & index1,string & index1Q,string & index2,string & index2Q,bool hasIndex1,bool hasIndex2);



public:
    BAMWriter(string filename);
    ~BAMWriter();
    void writeHeader();
    void writePairedSequence(string name,
			     string firstSeq,
			     string firstQual,
			     string secondSeq,
			     string secondQual,
			     string index1,
			     string index1Q,
			     string index2,
			     string index2Q,
			     bool hasIndex1,
			     bool hasIndex2);

    void writeSingleSequence(string name,
			     string firstSeq,
			     string firstQual,
			     string index1,
			     string index1Q,
			     string index2,
			     string index2Q,
			     bool hasIndex1,
			     bool hasIndex2);
    void close();
};
#endif
