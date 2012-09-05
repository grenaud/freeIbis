#ifndef _DBSEQ_H_
#define _DBSEQ_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include "param.h"

using namespace std;

struct OneBfa
{
	bit32_t n;  //count
	bit24_t *s;
};
struct RefTitle
{
	string name;
	bit32_t size;
};
struct Block
{
	bit32_t id;
	bit32_t begin;
	bit32_t end;
};
struct KmerLoc
{
	bit32_t n1, n2, n3; //ab, ac, ad seed
	ref_id_t *id1, *id2, *id3;
	ref_loc_t *loc1, *loc2, *loc3;
};

class RefSeq
{
public:
	RefSeq();
	ref_loc_t LoadNextSeq(ifstream &fin);
	void BinSeq(OneBfa &a);
	void UnmaskRegion();
	void Run_ConvertBinseq(ifstream &fin);
	inline bit32_t s_MakeSeed_1(bit24_t *_m, int _a);
	inline bit32_t s_MakeSeed_2(bit24_t *_m, int _a, int _b, int _c, int _d);
	
	void InitialIndex();
	void t_CalKmerFreq_ab();
	void t_CalKmerFreq_ac();
	void t_CalKmerFreq_ad();
	void AllocIndex();
	void t_CreateIndex_ab();
	void t_CreateIndex_ac();
	void t_CreateIndex_ad();
	void CreateIndex();
	void ReleaseIndex();
#ifdef THREAD
	void *t_CreateIndex_ab(void *);
	void *t_CreateIndex_ac(void *);
	void *t_CreateIndex_ad(void *);
#endif	
	
public:
	int total_num;
	bit64_t sum_length;
	vector<OneBfa> bfa;
	bit32_t total_kmers;
	KmerLoc *index;
	vector<RefTitle> title;	
protected:
	ref_id_t _count;
	string _name;
	string _seq;
	ref_loc_t _length;
public:	
	vector<Block> _blocks;  //unmasked ref region
};

#endif //_DBSEQ_H_
