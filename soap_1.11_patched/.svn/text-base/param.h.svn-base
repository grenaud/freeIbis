#ifndef _PARAM_H_
#define _PARAM_H_

#include<cstdlib>
#include<iostream>
#include<string>

using namespace std;

typedef unsigned char bit8_t;
typedef unsigned short bit16_t;
typedef unsigned bit32_t;
typedef unsigned long long bit64_t;
struct bit24_t
{
	unsigned a:24;
};

#ifdef DB_CHR  // seqs <256, length <4Gb
	typedef bit8_t ref_id_t;
	typedef bit32_t ref_loc_t;
#endif
#ifdef DB_CONTIG // seqs <65K, length <4Gb
	typedef bit16_t ref_id_t;
	typedef bit32_t ref_loc_t;
#endif
#ifdef DB_SHORT // seqs <4G, length <65K
	typedef bit32_t ref_id_t;
	typedef bit16_t ref_loc_t;
#endif
#ifdef DB_HUGE // seqs<4G, length <4G
	typedef bit32_t ref_id_t;
	typedef bit32_t ref_loc_t;
#endif

class Param
{
public:
	Param();
	void SetSeedSize(int n);
	void SetMrnaTag(int n);
	void BuildMismatchTable();
public:
	int num_procs;  //number of parallel processors
	
	int chains;   //0: both; 1: direct only; 2: complementary only
	//dbseq
	int max_dbseq_size;
	int append_dbseq_size;
	//read
	int read_size;
	int max_ns;     //throw out reads containning >=max_ns 'N's
	int trim_lowQ;  //trim low-quality at 3'-end, or not?
	//quality
	bit8_t zero_qual;
	bit8_t qual_threshold;
	bit8_t default_qual;
	//pair-end mapping
	int min_insert;
	int max_insert;
	int optimize_output_SV;  //if a pair cannot align with proper orientation and distance, very likely a strctural variation happen here. we prefer to report hit of read 'a' and 'b' with smallest distance, so that to help detect structural variations
	//seed
	int half_seed_size;
	int seed_size;
	bit32_t half_seed_bits;
	bit32_t seed_bits;
	int min_read_size;
	//alignment
	int max_snp_num;   //maximum number of snps on one read allowed
	int max_gap_size;  //maximum gap size in a read
	int gap_edge;   //will not allow gap exist at the boundaries of a read
	int max_num_hits;   //maximum number of equal best hits, smaller will be faster
	//report hits
	int report_repeat_hits;   //how report repeat hits? 0: no, 1: pick one randomly, 2: report all
	bool output_id;   //1: output read id, 1: out read index
	
	//for mRNA tag alignment
	int tag_type;
	string tag_seq;
	int tag_remain;
	
	//for miRNA alignment
	string adapter;
	int admis;
	int mirna_min;
	int mirna_max;
	//alphabet table
	string useful_nt;
	string nx_nt;
	//mismatch table
	bit8_t num_mismatch[0xffffff];
};

#endif //_PARAM_H_
