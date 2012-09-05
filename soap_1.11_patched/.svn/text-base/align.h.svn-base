#ifndef _ALIGN_H_
#define _ALIGN_H_

#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<string>
#include<algorithm>

#include "dbseq.h"
#include "param.h"
#include "reads.h"

using namespace std;

#ifdef READ_36
const int FIXSIZE=48; //read size <=36bp, then (3+1)*bit24_t.
const int FIXELEMENT=4;
#endif

#ifdef READ_48
const int FIXSIZE=60; //read size <=48bp, then (4+1)*bit24_t.
const int FIXELEMENT=5;
#endif

#ifdef READ_60
const int FIXSIZE=72; //read size <=60bp, then (5+1)*bit24_t.
const int FIXELEMENT=6;
#endif

#ifdef READ_72
const int FIXSIZE=84; //read size <=72bp, then (6+1)*bit24_t.
const int FIXELEMENT=7;
#endif

#ifdef READ_84
const int FIXSIZE=96; //read size <=96bp, then (7+1)*bit24_t.
const int FIXELEMENT=8;
#endif
 
const int MAXSNPS=5;

extern Param param;
extern char rev_char[];

struct SeedProfile
{
	bit8_t a;  //offset of part a on binary seq
	bit8_t b;  //part b
	bit8_t b1, b2;  //begin element when creating seed
	bit8_t s1, s2;  //shift when creating seed from binary seq
};

struct Hit
{
	ref_id_t chr;       //index of chr
	ref_loc_t loc;     //location of first bp on reference seq, count from 0
	bit8_t z;  /*No. of bseq, [0, 4);
						 if gapped hit, 100+n: insertion on query at site n;
						 200+n: deletion on query at site n;
						 */
};

/*
struct GapHit
{
	ref_id_t chr;
	ref_loc_t loc;
	bool del;  //0: insertion; 1: deletion, on query
	bit8_t offset;
};
*/
class SingleAlign
{
public:
	SingleAlign();
	void ImportFileFormat(int format);
	void ImportBatchReads(bit32_t n, vector<ReadInf> &a);
	int CountNs();
	int TrimLowQual();
	int TrimLowQual_2();
	int TrimLowQual_3();
	void ConvertBinaySeq();
	inline void GenerateSeeds_1(int n);
	inline void GenerateSeeds_2(int n);
	inline void GenerateSeeds_3(int n);
	inline int CountMismatch(bit24_t *q, bit24_t *r, bit24_t *s);
	inline bool UnequalTag_0(ref_id_t id, ref_loc_t loc, RefSeq &ref);
	inline bool UnequalTag_1(ref_id_t id, ref_loc_t loc, RefSeq &ref);
	void SnpAlign_0(RefSeq &ref);
	void SortExactHits(void);
	void SortExactcHits(void);
	void SnpAlign_1(RefSeq &ref);
	inline bool equal_loc(Hit a);
	void SnpAlign_2(RefSeq &ref);
	int SnpAlign_range(bool chain, ref_id_t id, ref_loc_t left_end, ref_loc_t right_end, RefSeq &ref);
	
	void DelBp(bit24_t *ori, bit24_t *ori_reg, bit32_t pos, bit32_t g);
	void InsBp(bit24_t *ori, bit24_t *ori_reg, bit32_t pos, bit32_t g);
	int GapAlign(RefSeq &ref);
	int GapAlign_range(bool chain, ref_id_t id, ref_loc_t left_end, ref_loc_t right_end, RefSeq &ref);
	void ClearHits();
	int RunAlign(RefSeq &ref);
	int FilterReads();
	int CountStringMismatch(int offset, string &s1, string s2);
	void Do_Batch(RefSeq &ref);
	void DetectSnpSites(bool chain, Hit *hit, RefSeq &ref);
	void SetFlag(char c);
	void StringAlign(RefSeq &ref, string &os);
	void Reverse_Seq();
	void Reverse_Qual();
	void s_OutHit(int chain, size_t n, bit8_t nspsn, Hit *hit, bool sig, RefSeq &ref, string &os);
	void s_OutGapHit(int chain, size_t n, bit8_t g, Hit *hit, RefSeq &ref, string &os);

public:
	int _format;
	//reads
	vector<ReadInf>::iterator _pread;
	string _ori_read_seq;
	string _ori_read_qual;
	string _revseq;
	string _revqual;
	bit32_t num_reads;
	vector<ReadInf> mreads;
	bit32_t n_aligned;
	//binary seq
	bit24_t bseq[12][FIXELEMENT];
	bit24_t reg[12][FIXELEMENT];
	bit24_t cbseq[12][FIXELEMENT];
	bit24_t creg[12][FIXELEMENT];
	
	//seed
	SeedProfile profile[6][4];
	bit32_t seeds[6][4];
	bit32_t cseeds[6][4];
	
	//for gapped align
	bit24_t gbseq[FIXELEMENT+1];
	bit24_t greg[FIXELEMENT+1];
	bit32_t _leftbits[13];
	bit32_t _rightbits[13];
	bit32_t _gc1[2][4];
	bit32_t _gc2[2][4];
	bit32_t _num1[2][4];
	bit32_t _num2[2][4];
	ref_id_t *_id1[2][4];
	ref_loc_t *_loc1[2][4];
	ref_id_t *_id2[2][4];
	ref_loc_t *_loc2[2][4];
	
	//for mRNA tag alignment
	bit24_t mrna_tag_seq[12][2];
	bit24_t mrna_tag_reg[12][2];
	bit24_t mrna_tag_cseq[12][2];
	
	//for miRNA alignment
	bit32_t mirna_adapter[16];
	bit32_t mirna_ada_reg[16];
	
	//alignment hits
	Hit hits[MAXSNPS+1][MAXHITS+1];
	Hit chits[MAXSNPS+1][MAXHITS+1];
	Hit bound_hits[MAXSNPS+1][MAXHITS+1];
	Hit gaphits[MAXHITS+1];
	Hit cgaphits[MAXHITS+1];
	Hit bound_gaphits[MAXHITS+1];
	
	int _cur_n_hit[MAXSNPS+1];
	int _cur_n_chit[MAXSNPS+1];
	int _cur_n_gaphit, _cur_n_cgaphit;
	int _cur_n_boundhit[MAXSNPS+1];
	int _cur_n_boundgaphit;
	
	bit8_t _gap_size;
	
	char _setname; //name of the set of data, like 'a'
	string _str_align;   //align results, prepare for output
protected:
	//local variables
	SeedProfile *_pro;
	bit32_t _seed;
	bit24_t _a, _b;
	string::iterator _sp;
	string::reverse_iterator _sq;
	Hit _hit;
	Hit _gaphit;
	string _str;
	ref_id_t *_refid;
	ref_loc_t *_refloc;
	ref_id_t *_lid, *_rid, *_tid;
	ref_loc_t *_lp, *_rp, *_tp;
	ref_loc_t _lb, _rb;
	
	int _tmp_n_hit[MAXSNPS+1];
	int _tmp_n_gaphit;
	int _tmp_n_cgaphit;
	
	char _ch[1000];
	vector<bit8_t> _sites;
};
	
/*n=0: ab; 1: cd; 2: bc; 3: ac; 4: bd; 5: ad*/
//n<3: ab, cd, bc	
inline void SingleAlign::GenerateSeeds_1(int n)
{
	_pro=profile[n];
	seeds[n][0]=_pro->s1>=24? (bseq[0][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (bseq[0][_pro->b1].a<<(24-_pro->s1)|bseq[0][_pro->b1+1].a>>_pro->s1)&param.seed_bits;
	cseeds[n][0]=_pro->s1>=24? (cbseq[0][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (cbseq[0][_pro->b1].a<<(24-_pro->s1)|cbseq[0][_pro->b1+1].a>>_pro->s1)&param.seed_bits;
	_pro++;
	seeds[n][1]=_pro->s1>=24? (bseq[1][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (bseq[1][_pro->b1].a<<(24-_pro->s1)|bseq[1][_pro->b1+1].a>>_pro->s1)&param.seed_bits;
	cseeds[n][1]=_pro->s1>=24? (cbseq[1][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (cbseq[1][_pro->b1].a<<(24-_pro->s1)|cbseq[1][_pro->b1+1].a>>_pro->s1)&param.seed_bits;
	_pro++;
	seeds[n][2]=_pro->s1>=24? (bseq[2][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (bseq[2][_pro->b1].a<<(24-_pro->s1)|bseq[2][_pro->b1+1].a>>_pro->s1)&param.seed_bits;
	cseeds[n][2]=_pro->s1>=24? (cbseq[2][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (cbseq[2][_pro->b1].a<<(24-_pro->s1)|cbseq[2][_pro->b1+1].a>>_pro->s1)&param.seed_bits;
	_pro++;
	seeds[n][3]=_pro->s1>=24? (bseq[3][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (bseq[3][_pro->b1].a<<(24-_pro->s1)|bseq[3][_pro->b1+1].a>>_pro->s1)&param.seed_bits;
	cseeds[n][3]=_pro->s1>=24? (cbseq[3][_pro->b1].a>>(_pro->s1-24))&param.seed_bits : (cbseq[3][_pro->b1].a<<(24-_pro->s1)|cbseq[3][_pro->b1+1].a>>_pro->s1)&param.seed_bits;

//	cout<<"build seed: 0  "<<seeds[n][0]<<endl;	
//	cout<<"build seed: 1  "<<seeds[n][1]<<endl;
//	cout<<"build seed: 2  "<<seeds[n][2]<<endl;
//	cout<<"build seed: 3  "<<seeds[n][3]<<endl;
}
//n=3 ac, 5 ad
inline void SingleAlign::GenerateSeeds_2(int n)
{	
	_pro=profile[n];
	seeds[n][0]=(((bseq[0][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? bseq[0][_pro->b2].a>>(_pro->s2-24): bseq[0][_pro->b2].a<<(24-_pro->s2)|bseq[0][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;
	cseeds[n][0]=(((cbseq[0][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? cbseq[0][_pro->b2].a>>(_pro->s2-24): cbseq[0][_pro->b2].a<<(24-_pro->s2)|cbseq[0][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;
	_pro++;
	seeds[n][1]=(((bseq[1][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? bseq[1][_pro->b2].a>>(_pro->s2-24): bseq[1][_pro->b2].a<<(24-_pro->s2)|bseq[1][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;
	cseeds[n][1]=(((cbseq[1][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? cbseq[1][_pro->b2].a>>(_pro->s2-24): cbseq[1][_pro->b2].a<<(24-_pro->s2)|cbseq[1][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;
	_pro++;
	seeds[n][2]=(((bseq[2][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? bseq[2][_pro->b2].a>>(_pro->s2-24): bseq[2][_pro->b2].a<<(24-_pro->s2)|bseq[2][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;
	cseeds[n][2]=(((cbseq[2][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? cbseq[2][_pro->b2].a>>(_pro->s2-24): cbseq[2][_pro->b2].a<<(24-_pro->s2)|cbseq[2][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;
	_pro++;
	seeds[n][3]=(((bseq[3][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? bseq[3][_pro->b2].a>>(_pro->s2-24): bseq[3][_pro->b2].a<<(24-_pro->s2)|bseq[3][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;
	cseeds[n][3]=(((cbseq[3][_pro->b1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size)
		|(_pro->s2>=24? cbseq[3][_pro->b2].a>>(_pro->s2-24): cbseq[3][_pro->b2].a<<(24-_pro->s2)|cbseq[3][_pro->b2+1].a>>_pro->s2) &param.half_seed_bits;																	
}
//n=4 bd
inline void SingleAlign::GenerateSeeds_3(int n)
{
	_pro=profile[n];
	seeds[n][0]=((_pro->s1>=24? bseq[0][_pro->b1].a>>(_pro->s1-24): bseq[0][_pro->b1].a<<(24-_pro->s1)|bseq[0][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? bseq[0][_pro->b2].a>>(_pro->s2-24): bseq[0][_pro->b2].a<<(24-_pro->s2)|bseq[0][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;
	cseeds[n][0]=((_pro->s1>=24? cbseq[0][_pro->b1].a>>(_pro->s1-24): cbseq[0][_pro->b1].a<<(24-_pro->s1)|cbseq[0][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? cbseq[0][_pro->b2].a>>(_pro->s2-24): cbseq[0][_pro->b2].a<<(24-_pro->s2)|cbseq[0][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;
	_pro++;
	seeds[n][1]=((_pro->s1>=24? bseq[1][_pro->b1].a>>(_pro->s1-24): bseq[1][_pro->b1].a<<(24-_pro->s1)|bseq[1][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? bseq[1][_pro->b2].a>>(_pro->s2-24): bseq[1][_pro->b2].a<<(24-_pro->s2)|bseq[1][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;
	cseeds[n][1]=((_pro->s1>=24? cbseq[1][_pro->b1].a>>(_pro->s1-24): cbseq[1][_pro->b1].a<<(24-_pro->s1)|cbseq[1][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? cbseq[1][_pro->b2].a>>(_pro->s2-24): cbseq[1][_pro->b2].a<<(24-_pro->s2)|cbseq[1][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;
	_pro++;
	seeds[n][2]=((_pro->s1>=24? bseq[2][_pro->b1].a>>(_pro->s1-24): bseq[2][_pro->b1].a<<(24-_pro->s1)|bseq[2][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? bseq[2][_pro->b2].a>>(_pro->s2-24): bseq[2][_pro->b2].a<<(24-_pro->s2)|bseq[2][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;
	cseeds[n][2]=((_pro->s1>=24? cbseq[2][_pro->b1].a>>(_pro->s1-24): cbseq[2][_pro->b1].a<<(24-_pro->s1)|cbseq[2][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? cbseq[2][_pro->b2].a>>(_pro->s2-24): cbseq[2][_pro->b2].a<<(24-_pro->s2)|cbseq[2][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;
	_pro++;
	seeds[n][3]=((_pro->s1>=24? bseq[3][_pro->b1].a>>(_pro->s1-24): bseq[3][_pro->b1].a<<(24-_pro->s1)|bseq[3][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? bseq[3][_pro->b2].a>>(_pro->s2-24): bseq[3][_pro->b2].a<<(24-_pro->s2)|bseq[3][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;
	cseeds[n][3]=((_pro->s1>=24? cbseq[3][_pro->b1].a>>(_pro->s1-24): cbseq[3][_pro->b1].a<<(24-_pro->s1)|cbseq[3][_pro->b1+1].a>>_pro->s1)&param.half_seed_bits)<<param.seed_size
		|(_pro->s2>=24? cbseq[3][_pro->b2].a>>(_pro->s2-24): cbseq[3][_pro->b2].a<<(24-_pro->s2)|cbseq[3][_pro->b2+1].a>>_pro->s2)&param.half_seed_bits;

}	

inline int SingleAlign::CountMismatch(bit24_t *q, bit24_t *r, bit24_t *s)
{
#ifdef READ_36
//	cout<<q->a<<"  "<<s->a<<"  "<<r->a<<endl;
//	cout<<(q+1)->a<<"  "<<(s+1)->a<<"  "<<(r+1)->a<<endl;
//	cout<<(q+2)->a<<"  "<<(s+2)->a<<"  "<<(r+2)->a<<endl;
//	cout<<(q+3)->a<<"  "<<(s+3)->a<<"  "<<(r+3)->a<<endl;
	return param.num_mismatch[(q->a^s->a)&r->a]+param.num_mismatch[((q+1)->a^(s+1)->a)&(r+1)->a]
		+param.num_mismatch[((q+2)->a^(s+2)->a)&(r+2)->a]+param.num_mismatch[((q+3)->a^(s+3)->a)&(r+3)->a];
#endif
#ifdef READ_48
	return param.num_mismatch[(q->a^s->a)&r->a]+param.num_mismatch[((q+1)->a^(s+1)->a)&(r+1)->a]
		+param.num_mismatch[((q+2)->a^(s+2)->a)&(r+2)->a]+param.num_mismatch[((q+3)->a^(s+3)->a)&(r+3)->a]
		+param.num_mismatch[((q+4)->a^(s+4)->a)&(r+4)->a];
#endif
#ifdef READ_60
	return param.num_mismatch[(q->a^s->a)&r->a]+param.num_mismatch[((q+1)->a^(s+1)->a)&(r+1)->a]
		+param.num_mismatch[((q+2)->a^(s+2)->a)&(r+2)->a]+param.num_mismatch[((q+3)->a^(s+3)->a)&(r+3)->a]
		+param.num_mismatch[((q+4)->a^(s+4)->a)&(r+4)->a]+param.num_mismatch[((q+5)->a^(s+5)->a)&(r+5)->a];
#endif
#ifdef READ_72
	return param.num_mismatch[(q->a^s->a)&r->a]+param.num_mismatch[((q+1)->a^(s+1)->a)&(r+1)->a]
		+param.num_mismatch[((q+2)->a^(s+2)->a)&(r+2)->a]+param.num_mismatch[((q+3)->a^(s+3)->a)&(r+3)->a]
		+param.num_mismatch[((q+4)->a^(s+4)->a)&(r+4)->a]+param.num_mismatch[((q+5)->a^(s+5)->a)&(r+5)->a]+param.num_mismatch[((q+6)->a^(s+6)->a)&(r+6)->a];
#endif
#ifdef READ_84
	return param.num_mismatch[(q->a^s->a)&r->a]+param.num_mismatch[((q+1)->a^(s+1)->a)&(r+1)->a]
		+param.num_mismatch[((q+2)->a^(s+2)->a)&(r+2)->a]+param.num_mismatch[((q+3)->a^(s+3)->a)&(r+3)->a]
		+param.num_mismatch[((q+4)->a^(s+4)->a)&(r+4)->a]+param.num_mismatch[((q+5)->a^(s+5)->a)&(r+5)->a]+param.num_mismatch[((q+6)->a^(s+6)->a)&(r+6)->a]+param.num_mismatch[((q+7)->a^(s+7)->a)&(r+7)->a];
#endif
}	
inline void SingleAlign::Reverse_Seq()
{
	_revseq=_pread->seq;
	reverse(_revseq.begin(), _revseq.end());
	for(string::iterator p=_revseq.begin(); p!=_revseq.end(); p++)
		*p=rev_char[*p];
}

inline void SingleAlign::Reverse_Qual()
{
	_revqual=_pread->qual;
	reverse(_revqual.begin(), _revqual.end());
}
#endif //_ALIGN_H_
