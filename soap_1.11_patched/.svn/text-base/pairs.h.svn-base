#ifndef _PAIRS_H_
#define _PAIRS_H_

#include<cmath>

#include "dbseq.h"
#include "reads.h"
#include "align.h"

using namespace std;

struct PairHit
{
	bool chain;
	bit8_t na, nb;   //# of snps
	Hit a;
	Hit b;
};

class PairAlign
{
public:
	PairAlign();
	void SetFlag();
	void ImportFileFormat(int format1, int format2);
	void ImportBatchReads(bit32_t n, vector<ReadInf> &a1, vector<ReadInf> &a2);
	int GetExactPairs();
	int GetExact2SnpPairs(RefSeq &ref);
	int GetSnp2SnpPairs(RefSeq &ref);
	int GetExact2GapPairs(RefSeq &ref);
	int RunAlign(RefSeq &ref);
	void Do_Batch(RefSeq &ref);
	void StringAlign(RefSeq &ref, string &os);
	void StringAlign_ClosestUnpair(RefSeq &ref, string &os);
	
public:	
	SingleAlign _sa;
	SingleAlign _sb;
	bit32_t num_reads;
	bit32_t n_aligned_pairs, n_aligned_a, n_aligned_b;	
	string _str_align;
	string _str_align_unpair;
protected:
	bit32_t _cur_n_hits[2*MAXSNPS+1];
	bit32_t _cur_n_gaphits[MAXGAP+1];
	PairHit pairhits[2*MAXSNPS+1][MAXHITS+1];
	PairHit pairgaphits[MAXGAP+1][MAXHITS+1];
};

#endif //_PAIR_ALIGH_H_
