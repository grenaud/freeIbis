#include "reads.h"

using namespace std;

extern Param param;

ReadClass::ReadClass()
{
	_index=0;
	mreads.resize(BatchNum);
}

void ReadClass::CheckFile(ifstream &fin)
{
	string s1,s2,s3,s4;
	char ch[1000];
	fin>>s1;
	fin.getline(ch, 1000);
	fin>>s2;
	fin.getline(ch, 1000);
	
	if('>' == s1[0])
		_file_format=1;
	else if('@' == s1[0]) {
		fin>>s3;		
		fin.getline(ch, 1000);
		fin>>s4;
		fin.getline(ch, 1000);
		_file_format=0;
		if(s2.size() != s4.size()) {
			cerr<<"fatal error: fq format, sequence length not equal to quality length\n";
			exit(1);
		}		
	}
	else {
		cerr<<"fatal error: unrecognizable format of reads file.\n";
		exit(1);
	}
	if(s2.size() < param.min_read_size) {
		cerr<<"fatal error: please set smaller seed size, <="<<2*((s2.size()-4+1)/4)<<endl;
		exit(1);
	}
	fin.seekg(0);
}

void ReadClass::InitialIndex()
{
	_index=0;
}
int ReadClass::LoadBatchReads(ifstream &fin)
{
	char ch[1000];
	char c;
	num=0;
	vector<ReadInf>::iterator p=mreads.begin();
	for(; num<BatchNum; p++,num++,_index++) {
		fin>>c;
		if(fin.eof())
			break;
		p->index=_index;
		fin>>p->name;
		fin.getline(ch,1000);
		fin>>p->seq;
		if(!_file_format) {//*.fq
			fin>>ch;
			fin.getline(ch, 1000);
			fin>>p->qual;
		}
		else
			p->qual=string(p->seq.size(), param.zero_qual+param.default_qual);
	}
	return num;
}

