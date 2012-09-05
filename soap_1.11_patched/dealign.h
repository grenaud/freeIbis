#ifndef _DEALIGN_H_
#define _DEALIGN_H_

#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<map>
#include<iomanip>

using namespace std;

typedef unsigned short bit16_t;
typedef unsigned int bit32_t;
typedef unsigned char bit8_t;
typedef unsigned long long bit64_t;

const int max_depth=65535;
const int DEFAULT_READ_LEN=32;
const int READID_LEN=32;
const int SEQ_LEN=84;
const int CHRID_LEN=32;
const int NOTE_LEN=50;

struct AlignInfo
{
	char id[READID_LEN];
	char seq[SEQ_LEN];
	char qual[SEQ_LEN];
	bit16_t nhits;
	char flag;
	bit16_t len;
	char chain;
	bit16_t chr;
	bit32_t loc;
	char note[NOTE_LEN];
};

class Dealign
{
public:
	void IniCount(char *ref, vector<vector<bit16_t> > &c);
	void OrderTitle();
	void CountHeadFreq(char *align_file, int flag);
	void CountDepth(char *align_file, int flag);
	void OutHist(char *out_file, vector<vector<bit16_t> > &c, int win, int slip, float zoom);
	void OutDistri(char *out_file, vector<vector<bit16_t> > &c, int win, int slip, float zoom);
	vector<bit8_t> GCcontent(string s, int win, int skip);
	void CalGC(char *ref, int win, int slip);
	void OutGCdot(char *out_file, vector<vector<bit16_t> > &c, int win, int slip, float zoom);
	void LoadChrOrder(char *chrorder_file);
	void LoadAlign(char *align_file);
	void SortAlign(int flag);
	void ExportAlign(char *align_file);
	void MergeSorted(char *infiles, char *outfile);
	void IniQC();
	void QC(char *align_file, int read_len, int flag, char zero_qual);
	void OutQC(char *out_file);

public:	
	vector<vector<bit16_t> > _head;
	vector<vector<bit16_t> > _depth;
protected:
	vector<string> _ref;
	vector<string> _reftitle;
	map<string, int> _titleorder;
	vector<AlignInfo> _align;
	vector<vector<bit8_t> > _gc;
	
	//QC
	bit32_t _ur[2][10];
	bit64_t _ntfreq[SEQ_LEN][5];  //A,C,G,T,N;   GC,AT
	bit64_t _mut[SEQ_LEN][6];     //A,C,G,T,# base,mismatch rate
	bit64_t _match_qual[SEQ_LEN][20];  //0~9,10~19,20~29,...
	bit64_t _mis_qual[SEQ_LEN][20];    //0~9,10~19,...
};

#endif //_DEALIGN_H_
