#include "pairs.h"

using namespace std;

extern Param param;

PairAlign::PairAlign()
{
	n_aligned_pairs=n_aligned_a=n_aligned_b=0;
	SetFlag();
}
void PairAlign::SetFlag()
{
	_sa.SetFlag('a');
	_sb.SetFlag('b');
}
void PairAlign::ImportFileFormat(int format1, int format2)
{
	_sa.ImportFileFormat(format1);
	_sb.ImportFileFormat(format2);
}
void PairAlign::ImportBatchReads(bit32_t n, vector<ReadInf> &a1, vector<ReadInf> &a2)
{
	_sa.ImportBatchReads(n, a1);
	_sb.ImportBatchReads(n, a2);
	num_reads=n;
}
int PairAlign::GetExactPairs()
{
	int i, j, k;
	PairHit pp;
	pp.na=pp.nb=0;
	//a+, b-
	pp.chain=0;
	if(_sa._cur_n_hit[0] && _sb._cur_n_chit[0]) {
		i=j=0;
		_sa.SortExactHits();
		_sb.SortExactcHits();
		while(i<_sa._cur_n_hit[0]) {
			while((j<_sb._cur_n_chit[0]) &&((_sb.chits[0][j].chr<_sa.hits[0][i].chr) 
				||((_sb.chits[0][j].chr==_sa.hits[0][i].chr) &&(_sb.chits[0][j].loc+_sb._pread->seq.size()<_sa.hits[0][i].loc+param.min_insert))))
				j++;
			k=j;
			while((k<_sb._cur_n_chit[0]) &&(_sb.chits[0][k].chr==_sa.hits[0][i].chr) 
				&&(_sb.chits[0][k].loc+_sb._pread->seq.size()<=_sa.hits[0][i].loc+param.max_insert)) {
				if(_cur_n_hits[0]>=param.max_num_hits)
					return 0;
				pp.a=_sa.hits[0][i];
				pp.b=_sb.chits[0][k];
				pairhits[0][_cur_n_hits[0]++]=pp;
				k++;
			}
			i++;
		}
	}
	//a-, b+
	pp.chain=1;
	if(_sa._cur_n_chit[0] && _sb._cur_n_hit[0]) {
		i=j=0;
		_sa.SortExactcHits();
		_sb.SortExactHits();		
		while(i<_sa._cur_n_chit[0]) {
			while((j<_sb._cur_n_hit[0]) &&((_sb.hits[0][j].chr<_sa.chits[0][i].chr) 
				||((_sb.hits[0][j].chr==_sa.chits[0][i].chr) &&(_sb.hits[0][j].loc<_sa.chits[0][i].loc+_sa._pread->seq.size()-param.max_insert))))
				j++;
			k=j;
			while((k<_sb._cur_n_hit[0]) &&(_sb.hits[0][k].chr==_sa.chits[0][i].chr) 
				&&(_sb.hits[0][k].loc<=_sa.chits[0][i].loc+_sa._pread->seq.size()-param.min_insert)) {
				if(_cur_n_hits[0]>=param.max_num_hits)
					return 0;				
				pp.a=_sa.chits[0][i];
				pp.b=_sb.hits[0][k];
				pairhits[0][_cur_n_hits[0]++]=pp;
				k++;
			}
			i++;
		}
	}
	if(_cur_n_hits[0]>0)
		return 0;
	return -1;	
}
int PairAlign::GetExact2SnpPairs(RefSeq &ref)
{
	PairHit pp;
	int i,j,nsnp;
	//a+ vs b-
	if(_sa._cur_n_hit[0]) {
		pp.chain=0;
		pp.na=0;
		for(i=0; i<_sa._cur_n_hit[0]; i++) {
			nsnp=_sb.SnpAlign_range(1, _sa.hits[0][i].chr, _sa.hits[0][i].loc+param.min_insert-_sb._pread->seq.size(), _sa.hits[0][i].loc+param.max_insert-_sb._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.a=_sa.hits[0][i];
				for(j=0; j<_sb._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=param.max_num_hits)
						break;
					pp.b=_sb.bound_hits[nsnp][j];
					pp.nb=nsnp;					
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}
		}
	}
	//a- vs b+
	if(_sa._cur_n_chit[0]) {
		pp.chain=1;
		pp.na=0;
		for(i=0; i<_sa._cur_n_chit[0]; i++) {
			nsnp=_sb.SnpAlign_range(0, _sa.chits[0][i].chr, (_sa.chits[0][i].loc<param.max_insert?0:(_sa.chits[0][i].loc-param.max_insert+_sa._pread->seq.size())), _sa.chits[0][i].loc-param.min_insert+_sa._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.a=_sa.chits[0][i];
				for(j=0; j<_sb._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=param.max_num_hits)
						break;		
					pp.b=_sb.bound_hits[nsnp][j];
					pp.nb=nsnp;			
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}			
		}	
	}
	//b+ vs a-
	if(_sb._cur_n_hit[0]) {
		pp.chain=1;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_hit[0]; i++) {
			nsnp=_sa.SnpAlign_range(1, _sb.hits[0][i].chr, _sb.hits[0][i].loc+param.min_insert-_sa._pread->seq.size(), _sb.hits[0][i].loc+param.max_insert-_sa._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.b=_sb.hits[0][i];
				for(j=0; j<_sa._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=param.max_num_hits)
						break;	
					pp.a=_sa.bound_hits[nsnp][j];
					pp.na=nsnp;				
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}			
		}		
	}
	//b- vs a+
	if(_sb._cur_n_chit[0]) {
		pp.chain=0;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_chit[0]; i++) {
			nsnp=_sa.SnpAlign_range(0, _sb.chits[0][i].chr, (_sb.chits[0][i].loc<param.max_insert?0:(_sb.chits[0][i].loc-param.max_insert+_sb._pread->seq.size())), _sb.chits[0][i].loc-param.min_insert+_sb._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.b=_sb.chits[0][i];
				for(j=0; j<_sa._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=param.max_num_hits)
						break;
					pp.a=_sa.bound_hits[nsnp][j];
					pp.na=nsnp;					
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}			
		}		
	}
	for(i=0; i<=param.max_snp_num; i++) {
		if(_cur_n_hits[i]>0)
			return i;
	}
	return -1;
}
int PairAlign::GetSnp2SnpPairs(RefSeq &ref)
{
	PairHit pp;
	int i,j,k,h,nsnp;
	for(i=1; i<=param.max_snp_num; i++) {
		//a+ vs b-
		if(_sa._cur_n_hit[i]) {
  		pp.chain=0;
  		pp.na=i;			
  		for(j=0; j<_sa._cur_n_hit[i]; j++) {
  			nsnp=_sb.SnpAlign_range(1, _sa.hits[i][j].chr, _sa.hits[i][j].loc+param.min_insert-_sb._pread->seq.size(), _sa.hits[i][j].loc+param.max_insert-_sb._pread->seq.size(), ref);
  			if(-1 !=nsnp) {
  				pp.a=_sa.hits[i][j];
  				for(k=0; k<_sb._cur_n_boundhit[nsnp]; k++) {
  					if(_cur_n_hits[nsnp+i]>=param.max_num_hits)
  						break;					
  					pp.b=_sb.bound_hits[nsnp][k];
  					pp.nb=nsnp;
  					pairhits[nsnp+i][_cur_n_hits[nsnp+i]++]=pp;
  				}
  			}		
  		}
  	}
		//a- vs b+
		if(_sa._cur_n_chit[i]) {
  		pp.chain=1;
  		pp.na=i;			
  		for(j=0; j<_sa._cur_n_chit[i]; j++) {
  			nsnp=_sb.SnpAlign_range(0, _sa.chits[i][j].chr, _sa.chits[i][j].loc-param.max_insert+_sa._pread->seq.size(), _sa.chits[i][j].loc-param.min_insert+_sa._pread->seq.size(), ref);
  			if(-1 !=nsnp) {
  				pp.a=_sa.chits[i][j];
  				for(k=0; k<_sb._cur_n_boundhit[nsnp]; k++) {
  					if(_cur_n_hits[nsnp+i]>=param.max_num_hits)
  						break;						
  					pp.b=_sb.bound_hits[nsnp][k];
  					pp.nb=nsnp;
  					pairhits[nsnp+i][_cur_n_hits[nsnp+i]++]=pp;
  				}
  			}			
  		}
  	}
		for(h=2; h<=i*2; h++) {
			if(_cur_n_hits[h]>0)
				return h;
		}	
	}
	return -1;
}
int PairAlign::GetExact2GapPairs(RefSeq &ref)
{
	PairHit pp;
	int i,j,k,h,gapsize;
	//a+ vs b-
	if(_sa._cur_n_hit[0]) {
		pp.chain=0;
		pp.na=0;
		for(i=0; i<_sa._cur_n_hit[0]; i++) {
			gapsize=_sb.GapAlign_range(1, _sa.hits[0][i].chr, _sa.hits[0][i].loc+param.min_insert-_sb._pread->seq.size(), _sa.hits[0][i].loc+param.max_insert-_sb._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.a=_sa.hits[0][i];
				for(j=0; j<_sb._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=param.max_num_hits)
						break;					
					pp.b=_sb.bound_gaphits[j];
					pp.nb=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}
	}
	//a- vs b+
	if(_sa._cur_n_chit[0]) {
		pp.chain=1;
		pp.na=0;
		for(i=0; i<_sa._cur_n_chit[0]; i++) {
			gapsize=_sb.GapAlign_range(0, _sa.chits[0][i].chr, _sa.chits[0][i].loc-param.max_insert+_sa._pread->seq.size(), _sa.chits[0][i].loc-param.min_insert+_sa._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.a=_sa.chits[0][i];
				for(j=0; j<_sb._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=param.max_num_hits)
						break;						
					pp.b=_sb.bound_gaphits[j];
					pp.nb=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}		
	}
	//b+ vs a-
	if(_sb._cur_n_hit[0]) {
		pp.chain=1;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_hit[0]; i++) {
			gapsize=_sa.GapAlign_range(1, _sb.hits[0][i].chr, _sb.hits[0][i].loc+param.min_insert-_sa._pread->seq.size(), _sb.hits[0][i].loc+param.max_insert-_sa._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.b=_sb.hits[0][i];
				for(j=0; j<_sa._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=param.max_num_hits)
						break;						
					pp.a=_sa.bound_gaphits[j];
					pp.na=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}
	}
	//b- vs a+
	if(_sb._cur_n_chit[0]) {
		pp.chain=0;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_chit[0]; i++) {
			gapsize=_sa.GapAlign_range(0, _sb.chits[0][i].chr, _sb.chits[0][i].loc-param.max_insert+_sb._pread->seq.size(), _sb.chits[0][i].loc-param.min_insert+_sb._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.b=_sb.chits[0][i];
				for(j=0; j<_sa._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=param.max_num_hits)
						break;						
					pp.a=_sa.bound_gaphits[j];
					pp.na=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}		
	}
	for(i=1; i<=MAXGAP; i++) {
		if(_cur_n_gaphits[i])
			return i;
	}
	return -1;
}
int PairAlign::RunAlign(RefSeq &ref)
{
	int i;
	for(i=0; i<=param.max_snp_num*2; i++)
		_cur_n_hits[i]=0;
	for(i=0; i<=MAXGAP; i++)
		_cur_n_gaphits[i]=0;
	_sa.ClearHits();
	_sb.ClearHits();
	_sa.ConvertBinaySeq();
	_sb.ConvertBinaySeq();
	//get exact+exact pairs
	_sa.GenerateSeeds_1(0);
	_sb.GenerateSeeds_1(0);
	_sa.SnpAlign_0(ref);
	_sb.SnpAlign_0(ref);
	if(((_sa._cur_n_hit[0]&&_sb._cur_n_chit[0]) ||(_sa._cur_n_chit[0]&&_sb._cur_n_hit[0])) &&(-1!=GetExactPairs()))
		return 1;
	//get exact+snp pairs
	if(param.max_snp_num>0) {
		_sa.GenerateSeeds_1(1);
		_sa.GenerateSeeds_1(2);
		_sa.GenerateSeeds_2(3);
		_sa.GenerateSeeds_2(5);
		_sa.GenerateSeeds_3(4);
		_sb.GenerateSeeds_1(1);
		_sb.GenerateSeeds_1(2);
		_sb.GenerateSeeds_2(3);
		_sb.GenerateSeeds_2(5);
		_sb.GenerateSeeds_3(4);
		if(-1!=GetExact2SnpPairs(ref))
			return 2;
		//snp alignment for a, then do snp align for b in the flanking region, get pairs
		_sa.SnpAlign_1(ref);
		_sa.SnpAlign_2(ref);
		if(-1!=GetSnp2SnpPairs(ref))
		{
			return 3;
		}	
	}
	//gap align in flanking region and get exact+gap pairs
	if(param.max_gap_size &&(_sa._cur_n_hit[0] ||_sb._cur_n_chit[0] ||_sa._cur_n_chit[0] ||_sb._cur_n_hit[0]) &&(-1!=GetExact2GapPairs(ref)))
	{
		return 4;
	}
	return 0;		
}
void PairAlign::Do_Batch(RefSeq &ref)
{
	bool detect_pairs;
	bool detect_a;
	bool detect_b;
	_str_align.clear();
	_str_align_unpair.clear();
	int tt=0;
	int filter1, filter2;
	if(0==param.trim_lowQ) {
		for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
			filter1=_sa.FilterReads();
			filter2=_sb.FilterReads();
			if(!(filter1 || filter2)) {
				if(RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned_pairs++;
				}
				else {
					detect_a=detect_b=0;
					if(_sa.RunAlign(ref)) {
						detect_a=1;
						n_aligned_a++;
					}
					if(_sb.RunAlign(ref)) {
						detect_b=1;
						n_aligned_b++;
					}
					if(detect_a && detect_b) {
						if(param.optimize_output_SV)
							StringAlign_ClosestUnpair(ref, _str_align_unpair);
						else {
							_sa.StringAlign(ref, _str_align_unpair);
							_sb.StringAlign(ref, _str_align_unpair);
						}
					}
					else if(detect_a)
						_sa.StringAlign(ref, _str_align_unpair);
					else if(detect_b)
						_sb.StringAlign(ref, _str_align_unpair);
				}
			}
			else {
				if(!filter1) {
					if(_sa.RunAlign(ref)) {
						_sa.StringAlign(ref, _str_align_unpair);
						n_aligned_a++;
					}
				}
				if(!filter2) {
					if(_sb.RunAlign(ref)) {
						_sb.StringAlign(ref, _str_align_unpair);
						n_aligned_b++;
					}
				}
			}
		}
	}
	else if(10>=param.trim_lowQ) {
		for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
			_sa._pread->seq.erase(_sa._pread->seq.size()-param.trim_lowQ, param.trim_lowQ);
			_sa._pread->qual.erase(_sa._pread->qual.size()-param.trim_lowQ, param.trim_lowQ);			
			_sb._pread->seq.erase(_sb._pread->seq.size()-param.trim_lowQ, param.trim_lowQ);		
			_sb._pread->qual.erase(_sb._pread->qual.size()-param.trim_lowQ, param.trim_lowQ);
			filter1=_sa.FilterReads();
			filter2=_sb.FilterReads();
			if(!(filter1 || filter2)) {
				if(RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned_pairs++;
				}
				else {
					detect_a=detect_b=0;
					if(_sa.RunAlign(ref)) {
						detect_a=1;
						n_aligned_a++;
					}
					if(_sb.RunAlign(ref)) {
						detect_b=1;
						n_aligned_b++;
					}
					if(detect_a && detect_b) {
						if(param.optimize_output_SV)
							StringAlign_ClosestUnpair(ref, _str_align_unpair);
						else {
							_sa.StringAlign(ref, _str_align_unpair);
							_sb.StringAlign(ref, _str_align_unpair);
						}
					}
					else if(detect_a)
						_sa.StringAlign(ref, _str_align_unpair);
					else if(detect_b)
						_sb.StringAlign(ref, _str_align_unpair);
				}
			}
			else {
				if(!filter1) {
					if(_sa.RunAlign(ref)) {
						_sa.StringAlign(ref, _str_align_unpair);
						n_aligned_a++;
					}
				}
				if(!filter2) {
					if(_sb.RunAlign(ref)) {
						_sb.StringAlign(ref, _str_align_unpair);
						n_aligned_b++;
					}
				}				
			}
		}		
	}
	else if(20>=param.trim_lowQ) {
		for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
			_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-10), param.trim_lowQ-10);
			_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-10), param.trim_lowQ-10);	
			_sa._pread->seq.erase(0,1);
			_sa._pread->qual.erase(0,1);		
			_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-10), param.trim_lowQ-10);		
			_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-10), param.trim_lowQ-10);
			_sb._pread->seq.erase(0,1);
			_sb._pread->qual.erase(0,1);
			
			filter1=_sa.FilterReads();
			filter2=_sb.FilterReads();
			if(!(filter1 || filter2)) {
				if(RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned_pairs++;
				}
				else {
					detect_a=detect_b=0;
					if(_sa.RunAlign(ref)) {
						detect_a=1;
						n_aligned_a++;
					}
					if(_sb.RunAlign(ref)) {
						detect_b=1;
						n_aligned_b++;
					}
					if(detect_a && detect_b) {
						if(param.optimize_output_SV)
							StringAlign_ClosestUnpair(ref, _str_align_unpair);
						else {
							_sa.StringAlign(ref, _str_align_unpair);
							_sb.StringAlign(ref, _str_align_unpair);
						}
					}
					else if(detect_a)
						_sa.StringAlign(ref, _str_align_unpair);
					else if(detect_b)
						_sb.StringAlign(ref, _str_align_unpair);
				}
			}
			else {
				if(!filter1) {
					if(_sa.RunAlign(ref)) {
						_sa.StringAlign(ref, _str_align_unpair);
						n_aligned_a++;
					}
				}
				if(!filter2) {
					if(_sb.RunAlign(ref)) {
						_sb.StringAlign(ref, _str_align_unpair);
						n_aligned_b++;
					}
				}				
			}
		}		
	}
	else if(30>=param.trim_lowQ) {
		for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
			//store the original seq and qual
			_sa._ori_read_seq = _sa._pread->seq;
			_sa._ori_read_qual = _sa._pread->qual;
			_sb._ori_read_seq = _sb._pread->seq;
			_sb._ori_read_qual = _sb._pread->qual;
			//try pair-end align first, if no hits, try single-read align later
			if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned_pairs++;
				continue;
			}
			_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
			_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
			_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
			_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
			if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned_pairs++;
				continue;
			}
			// if no hits, try single-read align later
			_sa._pread->seq = _sa._ori_read_seq;
			_sa._pread->qual = _sa._ori_read_qual;
			_sb._pread->seq = _sb._ori_read_seq;
			_sb._pread->qual = _sb._ori_read_qual;
			
			detect_a=detect_b=0;
			if(!_sa.FilterReads() && _sa.RunAlign(ref)) {
				detect_a=1;
				n_aligned_a++;
			}
			else {
				_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
				_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
				if(!_sa.FilterReads()) {
					if(_sa.RunAlign(ref)) {
						detect_a=1;
						n_aligned_a++;
					}
				}
			}
			if(!_sb.FilterReads() && _sb.RunAlign(ref)) {
				detect_b=1;
				n_aligned_b++;
			}
			else {
				_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
				_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
				if(!_sb.FilterReads()) {
					if(_sb.RunAlign(ref)) {
						detect_b=1;
						n_aligned_b++;
					}
				}
			}
			if(detect_a && detect_b) {
				if(param.optimize_output_SV)
					StringAlign_ClosestUnpair(ref, _str_align_unpair);
				else {
					_sa.StringAlign(ref, _str_align_unpair);
					_sb.StringAlign(ref, _str_align_unpair);
				}
			}
			else if(detect_a)
				_sa.StringAlign(ref, _str_align_unpair);
			else if(detect_b)
				_sb.StringAlign(ref, _str_align_unpair);			
		}		
	}
	else if(40>=param.trim_lowQ) {
		for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
			//store the original seq and qual
			_sa._ori_read_seq = _sa._pread->seq;
			_sa._ori_read_qual = _sa._pread->qual;
			_sb._ori_read_seq = _sb._pread->seq;
			_sb._ori_read_qual = _sb._pread->qual;
			//try pair-end align first, if no hits, try single-read align later
			if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned_pairs++;
				continue;
			}
			_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
			_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
			_sa._pread->seq.erase(0,1);
			_sa._pread->qual.erase(0,1);
			_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
			_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
			_sb._pread->seq.erase(0,1);
			_sb._pread->qual.erase(0,1);
			if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned_pairs++;
				continue;
			}
			// if no hits, try single-read align later
			_sa._pread->seq = _sa._ori_read_seq;
			_sa._pread->qual = _sa._ori_read_qual;
			_sb._pread->seq = _sb._ori_read_seq;
			_sb._pread->qual = _sb._ori_read_qual;
			
			detect_a=detect_b=0;
			if(!_sa.FilterReads() && _sa.RunAlign(ref)) {
				detect_a=1;
				n_aligned_a++;
			}
			else {
				_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
				_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
				_sa._pread->seq.erase(0,1);
				_sa._pread->qual.erase(0,1);				
				if(!_sa.FilterReads()) {
					if(_sa.RunAlign(ref)) {
						detect_a=1;
						n_aligned_a++;
					}
				}
			}
			if(!_sb.FilterReads() && _sb.RunAlign(ref)) {
				detect_b=1;
				n_aligned_b++;
			}
			else {
				_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
				_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
				_sb._pread->seq.erase(0,1);
				_sb._pread->qual.erase(0,1);				
				if(!_sb.FilterReads()) {
					if(_sb.RunAlign(ref)) {
						detect_b=1;
						n_aligned_b++;
					}
				}
			}
			if(detect_a && detect_b) {
				if(param.optimize_output_SV)
					StringAlign_ClosestUnpair(ref, _str_align_unpair);
				else {
					_sa.StringAlign(ref, _str_align_unpair);
					_sb.StringAlign(ref, _str_align_unpair);
				}
			}
			else if(detect_a)
				_sa.StringAlign(ref, _str_align_unpair);
			else if(detect_b)
				_sb.StringAlign(ref, _str_align_unpair);				
		}		
	}
	else if(50>=param.trim_lowQ) {
		for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
			//store the original seq and qual
			_sa._ori_read_seq = _sa._pread->seq;
			_sa._ori_read_qual = _sa._pread->qual;
			_sb._ori_read_seq = _sb._pread->seq;
			_sb._ori_read_qual = _sb._pread->qual;
			//try pair-end align first, if no hits, try single-read align later
			if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned_pairs++;
				continue;
			}
			detect_pairs=0;
			while(1) {
				_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
				_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
				_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
				_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
				if(_sa._pread->seq.size()<param.min_read_size)
					break;
				if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned_pairs++;
					detect_pairs=1;
					break;
				}
			}
			if(detect_pairs)
				continue;
			// if no hits, try single-read align later
			_sa._pread->seq = _sa._ori_read_seq;
			_sa._pread->qual = _sa._ori_read_qual;
			_sb._pread->seq = _sb._ori_read_seq;
			_sb._pread->qual = _sb._ori_read_qual;
			
			detect_a=detect_b=0;
			if(!_sa.FilterReads() && _sa.RunAlign(ref)) {
				_sa.StringAlign(ref, _str_align_unpair);
				detect_a=1;
				n_aligned_a++;
			}
			else {
				while(1) {
					_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
					_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
					if(_sa._pread->seq.size()<param.min_read_size)
						break;
					if(!_sa.FilterReads() && _sa.RunAlign(ref)) {
							n_aligned_a++;
							detect_a=1;
							break;
					}
				}
			}
			if(!_sb.FilterReads() && _sb.RunAlign(ref)) {
				detect_b=1;
				n_aligned_b++;
			}
			else {
				while(1) {
					_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
					_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
					if(_sb._pread->seq.size()<param.min_read_size)
						break;
					if(!_sb.FilterReads() && _sb.RunAlign(ref)) {
							n_aligned_b++;
							detect_b=1;
							break;
					}
				}
			}
			if(detect_a && detect_b) {
				if(param.optimize_output_SV)
					StringAlign_ClosestUnpair(ref, _str_align_unpair);
				else {
					_sa.StringAlign(ref, _str_align_unpair);
					_sb.StringAlign(ref, _str_align_unpair);
				}
			}
			else if(detect_a)
				_sa.StringAlign(ref, _str_align_unpair);
			else if(detect_b)
				_sb.StringAlign(ref, _str_align_unpair);				
		}		
	}
	else if(60>=param.trim_lowQ) {
		for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
			//store the original seq and qual
			_sa._ori_read_seq = _sa._pread->seq;
			_sa._ori_read_qual = _sa._pread->qual;
			_sb._ori_read_seq = _sb._pread->seq;
			_sb._ori_read_qual = _sb._pread->qual;
			//try pair-end align first, if no hits, try single-read align later
			if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned_pairs++;
				continue;
			}
			detect_pairs=0;
			_sa._pread->seq.erase(0,1);
			_sa._pread->qual.erase(0,1);
			_sb._pread->seq.erase(0,1);
			_sb._pread->qual.erase(0,1);
			while(1) {
				_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
				_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
				_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
				_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
				if(_sa._pread->seq.size()<param.min_read_size)
					break;
				if(!(_sa.FilterReads() || _sb.FilterReads()) && RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned_pairs++;
					detect_pairs=1;
					break;
				}
			}
			if(detect_pairs)
				continue;
			// if no hits, try single-read align later
			_sa._pread->seq = _sa._ori_read_seq;
			_sa._pread->qual = _sa._ori_read_qual;
			_sb._pread->seq = _sb._ori_read_seq;
			_sb._pread->qual = _sb._ori_read_qual;
			
			detect_a=detect_b=0;
			if(!_sa.FilterReads() && _sa.RunAlign(ref)) {
				detect_a=1;
				n_aligned_a++;
			}
			else {
				_sa._pread->seq.erase(0,1);
				_sa._pread->qual.erase(0,1);				
				while(1) {
					_sa._pread->seq.erase(_sa._pread->seq.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
					_sa._pread->qual.erase(_sa._pread->qual.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
					if(_sa._pread->seq.size()<param.min_read_size)
						break;
					if(!_sa.FilterReads() && _sa.RunAlign(ref)) {
							n_aligned_a++;
							detect_a=1;
							break;
					}
				}
			}
			if(!_sb.FilterReads() && _sb.RunAlign(ref)) {
				detect_b=1;
				n_aligned_b++;
			}
			else {
				_sb._pread->seq.erase(0,1);
				_sb._pread->qual.erase(0,1);
				while(1) {
					_sb._pread->seq.erase(_sb._pread->seq.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
					_sb._pread->qual.erase(_sb._pread->qual.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
					if(_sb._pread->seq.size()<param.min_read_size)
						break;
					if(!_sb.FilterReads() && _sb.RunAlign(ref)) {
							n_aligned_b++;
							detect_b=1;
							break;
					}
				}
			}
			if(detect_a && detect_b) {
				if(param.optimize_output_SV)
					StringAlign_ClosestUnpair(ref, _str_align_unpair);
				else {
					_sa.StringAlign(ref, _str_align_unpair);
					_sb.StringAlign(ref, _str_align_unpair);
				}
			}
			else if(detect_a)
				_sa.StringAlign(ref, _str_align_unpair);
			else if(detect_b)
				_sb.StringAlign(ref, _str_align_unpair);			
		}		
	}
//	cout<<_str_align<<endl;
}
void PairAlign::StringAlign(RefSeq &ref, string &os)
{
	_sa.Reverse_Seq();
	_sa.Reverse_Qual();
	_sb.Reverse_Seq();
	_sb.Reverse_Qual();	
	int i, j;
	//snp hits
	for(i=0; i<=param.max_snp_num*2; i++) {
		if(0==_cur_n_hits[i])
			continue;
		if(1==_cur_n_hits[i]) {
			if(pairhits[i][0].na==0)
				_sa.s_OutHit(pairhits[i][0].chain, 1, pairhits[i][0].na, &pairhits[i][0].a, 1, ref, os);
			else
				_sa.s_OutHit(pairhits[i][0].chain, 1, pairhits[i][0].na, &pairhits[i][0].a, 1, ref, os);
			if(pairhits[i][0].nb==0)
				_sb.s_OutHit(!pairhits[i][0].chain, 1, pairhits[i][0].nb, &pairhits[i][0].b, 1, ref, os);
			else
				_sb.s_OutHit(!pairhits[i][0].chain, 1, pairhits[i][0].nb, &pairhits[i][0].b, 1, ref, os);			
		}
		else if(1==param.report_repeat_hits) {   //randomly pick up one
			j=rand()%_cur_n_hits[i];
			_sa.s_OutHit(pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].na, &pairhits[i][j].a, 1, ref, os);
			_sb.s_OutHit(!pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].nb, &pairhits[i][j].b, 1, ref, os);
		}
		else if(2==param.report_repeat_hits) {   //output all repeat hits
			for(j=0; j<_cur_n_hits[i]; j++) {
				_sa.s_OutHit(pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].na, &pairhits[i][j].a, 1, ref, os);
				_sb.s_OutHit(!pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].nb, &pairhits[i][j].b, 1, ref, os);				
			}
		}
		return;
	}
	//gap hits
	for(i=1; i<=MAXGAP; i++) {
		if(0==_cur_n_gaphits[i])
			continue;
		if(1==_cur_n_gaphits[i]) {
			_sa.s_OutGapHit(pairgaphits[i][0].chain, _cur_n_gaphits[i], i, &pairgaphits[i][0].a, ref, os);
			_sb.s_OutGapHit(!pairgaphits[i][0].chain, _cur_n_gaphits[i], i, &pairgaphits[i][0].b, ref, os);
		}
		else if(1==param.report_repeat_hits) {
			j=rand()%_cur_n_gaphits[i];
			_sa.s_OutGapHit(pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].a, ref, os);
			_sb.s_OutGapHit(!pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].b, ref, os);
		}
		else if(2==param.report_repeat_hits) {
			for(j=0; j<_cur_n_gaphits[i]; j++) {						
				_sa.s_OutGapHit(pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].a, ref, os);
				_sb.s_OutGapHit(!pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].b, ref, os);
			}
		}
		return;
	}
}

void PairAlign::StringAlign_ClosestUnpair(RefSeq &ref, string &os)
{
	int nsnp_a=-1;
	int nsnp_b=-1;
	int chain_a, chain_b;
	int order_a, order_b;
	bit32_t score=0xffffffff;
	for(int i=0; i<=param.max_snp_num; i++) {
		if(_sa._cur_n_hit[i] || _sa._cur_n_chit[i]) {
			nsnp_a=i;
			break;
		}
	}
	for(int j=0; j<=param.max_snp_num; j++) {
		if(_sb._cur_n_hit[j] || _sb._cur_n_chit[j]) {
			nsnp_b=j;
			break;
		}
	}
	if((nsnp_a!=-1) && (nsnp_b!=-1)) {
		for(int i=0; i<_sa._cur_n_hit[nsnp_a]; i++) {
			chain_a=0;
			order_a=i;
			for(int j=0; j<_sb._cur_n_hit[nsnp_b]; j++) {
				chain_b=0;
				if((_sa.hits[nsnp_a][i].chr==_sb.hits[nsnp_b][j].chr)
					&& (abs(long(_sa.hits[nsnp_a][i].loc-_sb.hits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.hits[nsnp_a][i].loc-_sb.hits[nsnp_b][j].loc));
						order_b=j;
				}
			}
			for(int j=0; j<_sb._cur_n_chit[nsnp_b]; j++) {
				chain_b=1;
				if((_sa.hits[nsnp_a][i].chr==_sb.chits[nsnp_b][j].chr)
					&& (abs(long(_sa.hits[nsnp_a][i].loc-_sb.chits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.hits[nsnp_a][i].loc-_sb.chits[nsnp_b][j].loc));
						order_b=j;
				}
			}
		}
		for(int i=0; i<_sa._cur_n_chit[nsnp_a]; i++) {
			chain_a=1;
			order_a=i;
			for(int j=0; j<_sb._cur_n_hit[nsnp_b]; j++) {
				chain_b=0;
				if((_sa.chits[nsnp_a][i].chr==_sb.hits[nsnp_b][j].chr)
					&& (abs(long(_sa.chits[nsnp_a][i].loc-_sb.hits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.chits[nsnp_a][i].loc-_sb.hits[nsnp_b][j].loc));
						order_b=j;
				}
			}
			for(int j=0; j<_sb._cur_n_chit[nsnp_b]; j++) {
				chain_b=1;
				if((_sa.chits[nsnp_a][i].chr==_sb.chits[nsnp_b][j].chr)
					&& (abs(long(_sa.chits[nsnp_a][i].loc-_sb.chits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.chits[nsnp_a][i].loc-_sb.chits[nsnp_b][j].loc));
						order_b=j;
				}
			}			
		}
		if(score!=0xffffffff) {
			if(!chain_a)
				_sa.s_OutHit(0, _sa._cur_n_hit[nsnp_a]+_sa._cur_n_chit[nsnp_a], nsnp_a, &_sa.hits[nsnp_a][order_a], 1, ref, os);
			else
				_sa.s_OutHit(1, _sa._cur_n_hit[nsnp_a]+_sa._cur_n_chit[nsnp_a], nsnp_a, &_sa.chits[nsnp_a][order_a], 1, ref, os);
			if(!chain_b)
				_sb.s_OutHit(0, _sb._cur_n_hit[nsnp_b]+_sb._cur_n_chit[nsnp_b], nsnp_b, &_sb.hits[nsnp_b][order_b], 1, ref, os);
			else
				_sb.s_OutHit(1, _sb._cur_n_hit[nsnp_b]+_sb._cur_n_chit[nsnp_b], nsnp_b, &_sb.chits[nsnp_b][order_b], 1, ref, os);			
		}
		else {
			_sa.StringAlign(ref, os);
			_sb.StringAlign(ref, os);
		}
	}
	//a: gap; b: exact or snp
	else if((nsnp_a==-1) && (_sa._cur_n_gaphit || _sa._cur_n_cgaphit)) {
		for(int i=0; i<_sa._cur_n_gaphit; i++) {
			chain_a=0;
			order_a=i;
			for(int j=0; j<_sb._cur_n_hit[nsnp_b]; j++) {
				chain_b=0;
				if((_sa.gaphits[i].chr==_sb.hits[nsnp_b][j].chr)
					&& (abs(long(_sa.gaphits[i].loc-_sb.hits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.gaphits[i].loc-_sb.hits[nsnp_b][j].loc));
						order_b=j;
				}
			}
			for(int j=0; j<_sb._cur_n_chit[nsnp_b]; j++) {
				chain_b=1;
				if((_sa.gaphits[i].chr==_sb.chits[nsnp_b][j].chr)
					&& (abs(long(_sa.gaphits[i].loc-_sb.chits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.gaphits[i].loc-_sb.chits[nsnp_b][j].loc));
						order_b=j;
				}
			}
		}
		for(int i=0; i<_sa._cur_n_cgaphit; i++) {
			chain_a=0;
			order_a=i;
			for(int j=0; j<_sb._cur_n_hit[nsnp_b]; j++) {
				chain_b=0;
				if((_sa.cgaphits[i].chr==_sb.hits[nsnp_b][j].chr)
					&& (abs(long(_sa.cgaphits[i].loc-_sb.hits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.cgaphits[i].loc-_sb.hits[nsnp_b][j].loc));
						order_b=j;
				}
			}
			for(int j=0; j<_sb._cur_n_chit[nsnp_b]; j++) {
				chain_b=1;
				if((_sa.cgaphits[i].chr==_sb.chits[nsnp_b][j].chr)
					&& (abs(long(_sa.cgaphits[i].loc-_sb.chits[nsnp_b][j].loc))<score)) {
						score=abs(long(_sa.cgaphits[i].loc-_sb.chits[nsnp_b][j].loc));
						order_b=j;
				}
			}			
		}
		if(score!=0xffffffff) {
			if(!chain_a)
				_sa.s_OutGapHit(0, _sa._cur_n_gaphit+_sa._cur_n_cgaphit, _sa._gap_size, &_sa.gaphits[order_a], ref, os);
			else
				_sa.s_OutGapHit(1, _sa._cur_n_gaphit+_sa._cur_n_cgaphit, _sa._gap_size, &_sa.cgaphits[order_a], ref, os);
			if(!chain_b)
				_sb.s_OutHit(0, _sb._cur_n_hit[nsnp_b]+_sb._cur_n_chit[nsnp_b], nsnp_b, &_sb.hits[nsnp_b][order_b], 1, ref, os);
			else
				_sb.s_OutHit(1, _sb._cur_n_hit[nsnp_b]+_sb._cur_n_chit[nsnp_b], nsnp_b, &_sb.chits[nsnp_b][order_b], 1, ref, os);
		}	
	}
	//a: exact or snp; b: gap
	else if((nsnp_b==-1) && (_sb._cur_n_gaphit || _sb._cur_n_cgaphit)) {
		for(int i=0; i<_sb._cur_n_gaphit; i++) {
			chain_b=0;
			order_b=i;
			for(int j=0; j<_sa._cur_n_hit[nsnp_a]; j++) {
				chain_a=0;
				if((_sb.gaphits[i].chr==_sa.hits[nsnp_a][j].chr)
					&& (abs(long(_sb.gaphits[i].loc-_sa.hits[nsnp_a][j].loc))<score)) {
						score=abs(long(_sb.gaphits[i].loc-_sa.hits[nsnp_a][j].loc));
						order_a=j;
				}
			}
			for(int j=0; j<_sa._cur_n_chit[nsnp_a]; j++) {
				chain_a=1;
				if((_sb.gaphits[i].chr==_sa.chits[nsnp_a][j].chr)
					&& (abs(long(_sb.gaphits[i].loc-_sa.chits[nsnp_a][j].loc))<score)) {
						score=abs(long(_sb.gaphits[i].loc-_sa.chits[nsnp_a][j].loc));
						order_a=j;
				}
			}
		}
		for(int i=0; i<_sb._cur_n_cgaphit; i++) {
			chain_b=1;
			order_b=i;
			for(int j=0; j<_sa._cur_n_hit[nsnp_a]; j++) {
				chain_a=0;
				if((_sb.cgaphits[i].chr==_sa.hits[nsnp_a][j].chr)
					&& (abs(long(_sb.cgaphits[i].loc-_sa.hits[nsnp_a][j].loc))<score)) {
						score=abs(long(_sb.cgaphits[i].loc-_sa.hits[nsnp_a][j].loc));
						order_a=j;
				}
			}
			for(int j=0; j<_sa._cur_n_chit[nsnp_a]; j++) {
				chain_a=1;
				if((_sb.cgaphits[i].chr==_sa.chits[nsnp_a][j].chr)
					&& (abs(long(_sb.cgaphits[i].loc-_sa.chits[nsnp_a][j].loc))<score)) {
						score=abs(long(_sb.cgaphits[i].loc-_sa.chits[nsnp_a][j].loc));
						order_a=j;
				}
			}
		}
		if(score!=0xffffffff) {
			if(!chain_b)
				_sb.s_OutGapHit(0, _sb._cur_n_gaphit+_sb._cur_n_cgaphit, _sb._gap_size, &_sb.gaphits[order_b], ref, os);
			else
				_sb.s_OutGapHit(1, _sb._cur_n_gaphit+_sb._cur_n_cgaphit, _sb._gap_size, &_sb.cgaphits[order_b], ref, os);
			if(!chain_a)
				_sa.s_OutHit(0, _sa._cur_n_hit[nsnp_a]+_sa._cur_n_chit[nsnp_a], nsnp_a, &_sa.hits[nsnp_a][order_a], 1, ref, os);
			else
				_sa.s_OutHit(1, _sa._cur_n_hit[nsnp_a]+_sa._cur_n_chit[nsnp_a], nsnp_a, &_sa.chits[nsnp_a][order_a], 1, ref, os);
		}	
	}
	return;
}
