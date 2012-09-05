#ifndef _READS_H_
#define _READS_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include "param.h"

using namespace std;

const int BatchNum=10000;

struct ReadInf
{
	bit32_t index;
	string name;
	string seq;
	string qual;
};

class ReadClass
{
public:
	ReadClass();
	void CheckFile(ifstream &fin);
	void InitialIndex();
	int LoadBatchReads(ifstream &fin);
public:
	vector<ReadInf> mreads;
	bit32_t num;
	
	int _file_format;  //0: fq; 1: fa;
	bit32_t _index;
};

#endif //_READS_H_
