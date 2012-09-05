#include "align.h"

using namespace std;

extern bit8_t alphabet[];
extern bit8_t reg_alphabet[];
extern bit8_t rev_alphabet[];
extern char rev_char[];
extern char chain_flag[];
extern char nt_code[];
extern char revnt_code[];

//create seed profile
SingleAlign::SingleAlign()
{
	//
	if(param.max_snp_num >MAXSNPS) {
		cerr<<"fatal error, set smaller max_snp_num. (<=MAXSNPS)\n";
		exit(1);
	}
	n_aligned=0;
	//make seed profile
	/*n=0: ab; 1: cd; 2: bc; 3: ac; 4: bd; 5: ad*/
	/*
	      ==============
	    0 |---                  
	    1    |---
	    2       |---
	    3          |---
	*/
	SeedProfile e;
	for(int i=0; i<4; i++) {  //for 4 binary seqs
		//ab
		e.a=(i+3)&0xFC;  //start bp of seed part_a and part_b
		e.b1=e.a/12;    //start element of seed part_a and part_b
		e.s1=48-(param.seed_size+e.a%12)*2; //shift s1 to right for 2 elements
		profile[0][i]=e;
		//cd
		e.a=(param.half_seed_size*2+i+3)&0xFC; //[half_seed_size*2+i, half_seed_size*2+i+4), 4n
		e.b1=e.a/12;
		e.s1=48-(param.seed_size+e.a%12)*2;		//shift s1 to right for 2 elements
		profile[1][i]=e;
		//bc
		e.a=(param.half_seed_size+i+3)&0xFC;    //[half_seed_size+i, half_seed_size+i+4), 4n
		e.b1=e.a/12;
		e.s1=48-(param.seed_size+e.a%12)*2;   //for 2 elements
		profile[2][i]=e;				
		//ac
		e.a=(i+3)&0xFC; e.b=e.a+param.seed_size;
		e.b1=e.a/12; e.b2=e.b/12;
		e.s1=24-(param.half_seed_size+e.a%12)*2; //shift s1 to right for 1 element
		e.s2=48-(param.half_seed_size+e.b%12)*2;		
		profile[3][i]=e;
		//bd
		e.a=(param.half_seed_size+i+3)&0xFC; e.b=e.a+param.seed_size;
		e.b1=e.a/12; e.b2=e.b/12;
		e.s1=48-(param.half_seed_size+e.a%12)*2;
		e.s2=48-(param.half_seed_size+e.b%12)*2;		
		profile[4][i]=e;
		//ad
		e.a=(i+3)&0xFC; e.b=e.a+param.half_seed_size*3;
		e.b1=e.a/12; e.b2=e.b/12;
		e.s1=24-(param.half_seed_size+e.a%12)*2;  //for 1 element
		e.s2=48-(param.half_seed_size+e.b%12)*2;			
		profile[5][i]=e;
	}
	//make bits
	for(bit32_t i=0; i<=12; i++) {
		_rightbits[i]=(1<<(2*i))-1;
		_leftbits[i]=_rightbits[i]<<(24-2*i);
	}
	//make bits for mRNA alignment
	if(param.tag_type!=-1) {
		mrna_tag_seq[0][0].a=mrna_tag_seq[0][1].a=mrna_tag_reg[0][0].a=mrna_tag_reg[0][1].a=0;
		mrna_tag_cseq[0][0].a=mrna_tag_cseq[0][1].a=0;
		int i,j;
		for(i=0; i<param.tag_seq.size(); i++) {
			mrna_tag_seq[0][0].a |= alphabet[param.tag_seq[i]]<<(22-2*i);
			mrna_tag_reg[0][0].a |= 0x3<<(24-2*i);
		}
		for(i=param.tag_seq.size()-1, j=0; i>=0; i--, j++) {
			mrna_tag_cseq[0][0].a |= rev_alphabet[param.tag_seq[i]]<<(22-2*j);
		}
		for(i=1; i<12; i++) {
			mrna_tag_seq[i][0].a = mrna_tag_seq[i-1][0].a>>2;
			mrna_tag_cseq[i][0].a = mrna_tag_cseq[i-1][0].a>>2;
			mrna_tag_reg[i][0].a = mrna_tag_reg[i-1][0].a>>2;
			mrna_tag_seq[i][1].a = (mrna_tag_seq[i-1][1].a>>2)|((mrna_tag_seq[i-1][0].a&0x3)<<22);
			mrna_tag_cseq[i][1].a = (mrna_tag_cseq[i-1][1].a>>2)|((mrna_tag_cseq[i-1][0].a&0x3)<<22);
			mrna_tag_reg[i][1].a = (mrna_tag_reg[i-1][1].a>>2)|((mrna_tag_reg[i-1][0].a&0x3)<<22);	
		}
	}
	//make bits for miRNA alignment
	if(param.adapter.size()) {
		mirna_adapter[0]=0;
		mirna_ada_reg[0]=0;
		for(int i=0; (i<16)&&(i<param.adapter.size()); i++) {
			mirna_adapter[0] |= alphabet[param.adapter[i]]<<(32-2*i);
			mirna_ada_reg[0] |= 0x3<<(32-2*i);
		}
		for(int i=1; i<16; i++) {
			mirna_adapter[i] = mirna_adapter[i-1]>>2;
			mirna_ada_reg[i] = mirna_ada_reg[i-1]>>2;
		}
	}
}
	
void SingleAlign::ImportFileFormat(int format)
{
	_format=format;
}

void SingleAlign::ImportBatchReads(bit32_t n, vector<ReadInf> &a)
{
	num_reads=n;
	mreads=a;
}

int SingleAlign::CountNs()
{
	int n=0;
	for(_sp=_pread->seq.begin(); _sp!=_pread->seq.end(); _sp++)
		if(!reg_alphabet[*_sp])
			n++;
	return n;	
}

//trim low quality at 3'-end, cut at 3bp continuous high-quality bps
int SingleAlign::TrimLowQual()
{
	int n=0;
	int i=_pread->qual.size()+3;
	_sq=_pread->qual.rbegin();
	for(; _sq!=_pread->qual.rend(); i--,_sq++) {
		if(*_sq >=param.zero_qual+param.qual_threshold) {
			n++;
			if(n >= 3) { //define boundary at 3'-end of >=3 continuous high-quality bps
				if(i >= param.min_read_size) {
					_pread->qual.resize(i);
					_pread->seq.resize(i);
					return 1;
				}
			}
		}
		else
			n=0;
	}
	return 0;
}
//trim low quality, start from 5'-end, stop at accumulating 2 low-quality bps
int SingleAlign::TrimLowQual_2()
{
	int n=0;
	int i=0;
	for(_sp=_pread->qual.begin(); _sp!=_pread->qual.end(); _sp++, i++) {
		if(*_sp<param.zero_qual+param.qual_threshold) {
			n++;
			if(n>1)
				if(i>=param.min_read_size) {
					_pread->qual.resize(i);
					_pread->seq.resize(i);
				}
			else
				return 0;
		}
	}
	return 1;
}
//trim low quality at both ends, to get continuous high-quality fragments
int SingleAlign::TrimLowQual_3()
{
	return 1;
}
//shift 1bp to left
inline void RightShiftBinSeq(bit24_t *orib, bit24_t *newb)
{
#ifdef READ_36
	newb[0].a=(orib[0].a>>2);
	newb[1].a=(orib[1].a>>2)|((orib[0].a&0x3)<<22);
	newb[2].a=(orib[2].a>>2)|((orib[1].a&0x3)<<22);
	newb[3].a=(orib[3].a>>2)|((orib[2].a&0x3)<<22);
#endif	
#ifdef READ_48
	newb[0].a=(orib[0].a>>2);
	newb[1].a=(orib[1].a>>2)|((orib[0].a&0x3)<<22);
	newb[2].a=(orib[2].a>>2)|((orib[1].a&0x3)<<22);
	newb[3].a=(orib[3].a>>2)|((orib[2].a&0x3)<<22);
	newb[4].a=(orib[4].a>>2)|((orib[3].a&0x3)<<22);
#endif
#ifdef READ_60
	newb[0].a=(orib[0].a>>2);
	newb[1].a=(orib[1].a>>2)|((orib[0].a&0x3)<<22);
	newb[2].a=(orib[2].a>>2)|((orib[1].a&0x3)<<22);
	newb[3].a=(orib[3].a>>2)|((orib[2].a&0x3)<<22);
	newb[4].a=(orib[4].a>>2)|((orib[3].a&0x3)<<22);
	newb[5].a=(orib[5].a>>2)|((orib[4].a&0x3)<<22);
#endif
#ifdef READ_72
        newb[0].a=(orib[0].a>>2);
        newb[1].a=(orib[1].a>>2)|((orib[0].a&0x3)<<22);
        newb[2].a=(orib[2].a>>2)|((orib[1].a&0x3)<<22);
        newb[3].a=(orib[3].a>>2)|((orib[2].a&0x3)<<22);
        newb[4].a=(orib[4].a>>2)|((orib[3].a&0x3)<<22);
        newb[5].a=(orib[5].a>>2)|((orib[4].a&0x3)<<22);
        newb[6].a=(orib[6].a>>2)|((orib[5].a&0x3)<<22);
#endif
#ifdef READ_84
        newb[0].a=(orib[0].a>>2);
        newb[1].a=(orib[1].a>>2)|((orib[0].a&0x3)<<22);
        newb[2].a=(orib[2].a>>2)|((orib[1].a&0x3)<<22);
        newb[3].a=(orib[3].a>>2)|((orib[2].a&0x3)<<22);
        newb[4].a=(orib[4].a>>2)|((orib[3].a&0x3)<<22);
        newb[5].a=(orib[5].a>>2)|((orib[4].a&0x3)<<22);
        newb[6].a=(orib[6].a>>2)|((orib[5].a&0x3)<<22);
        newb[7].a=(orib[7].a>>2)|((orib[6].a&0x3)<<22);
#endif
}

//convert string seq to binary type
void SingleAlign::ConvertBinaySeq()
{	
	int i,j,h,g;
	//direct chain
	h=_a.a=_b.a=0;
	for(_sp=_pread->seq.begin(),i=1; _sp!=_pread->seq.end(); _sp++,i++) {
		_a.a<<=2; _b.a<<=2;
		_a.a|=alphabet[*_sp];
		_b.a|=reg_alphabet[*_sp];		
		if(0==i%12) {
			bseq[0][h]=_a;
			reg[0][h++]=_b;
			_a.a=_b.a=0;
		}		
	}
	for(; i!=FIXSIZE+1; i++) {
		_a.a<<=2; _b.a<<=2;
		if(0==i%12) {
			bseq[0][h]=_a;
			reg[0][h++]=_b;
			_a.a=_b.a=0;
		}		
	}
		i=0;
//		cout<<"bin seq: "<<i<<"  "<<bseq[i][0].a<<" "<<bseq[i][1].a<<" "<<bseq[i][2].a<<" "<<bseq[i][3].a<<endl;
//		cout<<"bin reg: "<<i<<"  "<<reg[i][0].a<<" "<<reg[i][1].a<<" "<<reg[i][2].a<<" "<<reg[i][3].a<<endl;	
	for(i=1; i!=12; i++) {
		RightShiftBinSeq(bseq[i-1], bseq[i]);
		RightShiftBinSeq(reg[i-1], reg[i]);
//		cout<<"bin seq: "<<i<<"  "<<bseq[i][0].a<<" "<<bseq[i][1].a<<" "<<bseq[i][2].a<<" "<<bseq[i][3].a<<endl;
//		cout<<"bin reg: "<<i<<"  "<<reg[i][0].a<<" "<<reg[i][1].a<<" "<<reg[i][2].a<<" "<<reg[i][3].a<<endl;
	}

	//reverse seq
	h=_a.a=_b.a=0;
	for(_sq=_pread->seq.rbegin(),i=1; _sq!=_pread->seq.rend(); _sq++,i++) {
		_a.a<<=2; _b.a<<=2;
		_a.a|=rev_alphabet[*_sq];
		_b.a|=reg_alphabet[*_sq];
		if(0==i%12) {
			cbseq[0][h]=_a;
			creg[0][h++]=_b;
			_a.a=_b.a=0;
		}		
	}
	for(; i!=FIXSIZE+1; i++) {
		_a.a<<=2; _b.a<<=2;
		if(0==i%12) {
			cbseq[0][h]=_a;
			creg[0][h++]=_b;
			_a.a=_b.a=0;
		}		
	}
	
	for(i=1; i!=12; i++) {
		RightShiftBinSeq(cbseq[i-1], cbseq[i]);
		RightShiftBinSeq(creg[i-1], creg[i]);
	}		
}

bool HitComp(Hit a, Hit b) 
{
	if(a.chr<b.chr)
		return 1;
	else if(a.chr==b.chr) {
		if(a.loc<b.loc)
			return 1;
		else
			return 0;
	}
	else
		return 0;
};

bool HitExist(Hit &a, Hit *pl, Hit *pr)
{
	Hit *pmiddle;
	while(pr>pl+1) {
		if((pl->chr==a.chr)&&(pl->loc==a.loc))
			return 1;
		pmiddle=pl+(pr-pl)/2;
		if((pmiddle->chr <a.chr)||((pmiddle->chr==a.chr) &&(pmiddle->loc<=a.loc)))
			pl=pmiddle;
		else
			pr=pmiddle;
	}
	if((pl->chr==a.chr)&&(pl->loc==a.loc))
		return 1;
	return 0;
};

//direct chain
inline bool SingleAlign::UnequalTag_0(ref_id_t id, ref_loc_t loc, RefSeq &ref)
{
	bit32_t m=(loc-param.tag_seq.size())/12;
	bit32_t n=(loc-param.tag_seq.size())%12;
	if(((ref.bfa[id].s[m].a^mrna_tag_seq[n][0].a)&mrna_tag_reg[n][0].a) 
		||((ref.bfa[id].s[m+1].a^mrna_tag_seq[n][1].a)&mrna_tag_reg[n][1].a))
		return 1;
	return 0;
}
//complementary chain
inline bool SingleAlign::UnequalTag_1(ref_id_t id, ref_loc_t loc, RefSeq &ref)
{
	bit32_t m=(loc+_pread->seq.size())/12;
	bit32_t n=(loc+_pread->seq.size())%12;
	if(((ref.bfa[id].s[m].a^mrna_tag_cseq[n][0].a)&mrna_tag_reg[n][0].a) 
		||((ref.bfa[id].s[m+1].a^mrna_tag_cseq[n][1].a)&mrna_tag_reg[n][1].a))
		return 1;
	return 0;	
}
/*n=0: ab; 1: bc; 2: cd; 3: ac; 4: bd; 5: ad*/
//for seed ab
void SingleAlign::SnpAlign_0(RefSeq &ref)
{
	bit32_t i,j,g,h,m,w,d;
	//direct chain
	if((0==param.chains) ||(1==param.chains)) {
	for(i=0; i!=4; i++) {
		_seed=seeds[0][i];
//		cout<<"direct: "<<i<<"  seed:  "<<_seed<<"  "<<ref.index[_seed].n1<<endl;
		if((m=ref.index[_seed].n1) ==0) continue;  //no match
		_refid=ref.index[_seed].id1;
		_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
		h=profile[0][i].b1;		 //# of units at left of seed on bseq
		for(j=0; j!=m; j++) {
			_hit.chr=_refid[j];
			_hit.loc=(_refloc[j]<<2)-profile[0][i].a+i;
//			cout<<"hit:  "<<_refloc[j]<<"   "<<j<<endl;
			if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
				continue;
			//special for mRNA tag alignment
			if((param.tag_type!=-1) && ((_hit.loc<4) || UnequalTag_0(_hit.chr, _hit.loc, ref)))
				continue;
			d=_hit.loc%12;
			w=CountMismatch(bseq[d], reg[d], ref.bfa[_hit.chr].s+_hit.loc/12);
//			cout<<"hit:  "<<_hit.loc<<"  "<<w<<endl;
			if(w >param.max_snp_num)
				continue;
			if(_cur_n_hit[w]>=param.max_num_hits) {
				if(w==0)
					break;
				else
					continue;
			}
			_hit.z=d;
			hits[w][_cur_n_hit[w]++]=_hit;
//			cout<<"+  "<<(int)_hit.chr<<"  "<<_hit.loc<<endl;
		}
	}
	}
	//complementary chain
	if((0==param.chains) ||(2==param.chains)) {
	for(i=0; i!=4; i++) {
		_seed=cseeds[0][i];
//		cout<<"compl: "<<i<<"  seed:  "<<_seed<<"  "<<ref.index[_seed].n1<<endl;
		if((m=ref.index[_seed].n1) ==0) continue;  //no match
//		cout<<i<<"  seed:  "<<_seed<<"  "<<ref.index[_seed].n1<<endl;
		_refid=ref.index[_seed].id1;
		_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
		h=profile[0][i].b1;		 //# of bytes at left of seed on bseq
		for(j=0; j!=m; j++) {
			_hit.chr=_refid[j];
			_hit.loc=(_refloc[j]<<2)-profile[0][i].a+i;
			if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
				continue;	
			//special for mRNA tag alignment
			if((param.tag_type!=-1) && UnequalTag_1(_hit.chr, _hit.loc, ref))
				continue;	
			d=_hit.loc%12;		
			w=CountMismatch(cbseq[d], creg[d], ref.bfa[_hit.chr].s+_hit.loc/12);
//			cout<<"hit:  "<<_hit.loc<<"  "<<w<<endl;			
			if(w >param.max_snp_num)
				continue;
			if(_cur_n_chit[w]>=param.max_num_hits) {
				if(w==0)
					break;
				else
					continue;
			}					
			_hit.z=d;
			chits[w][_cur_n_chit[w]++]=_hit;
		}
	}
	}
}
void SingleAlign::SortExactHits(void)
{
	sort(hits[0], hits[0]+_cur_n_hit[0], HitComp);
}
void SingleAlign::SortExactcHits(void)
{
	sort(chits[0], chits[0]+_cur_n_chit[0], HitComp);
}
//for seed cd
void SingleAlign::SnpAlign_1(RefSeq &ref)
{
	bit32_t i,j,g,h,m,w,d;	
	//direct chain
	if((0==param.chains) ||(1==param.chains)) {
	for(i=0; i!=param.max_snp_num+1; i++) {
		_tmp_n_hit[i]=_cur_n_hit[i];
		if(_cur_n_hit[i])
			sort(hits[i], hits[i]+_cur_n_hit[i], HitComp);
	}	
	for(i=0; i!=4; i++) {
		_seed=seeds[1][i];
//		cout<<"direct: "<<i<<"  seed:  "<<_seed<<"  "<<ref.index[_seed].n1<<endl;
		if((m=ref.index[_seed].n1) ==0) continue;  //no match
		_refid=ref.index[_seed].id1;
		_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
		h=profile[1][i].b1;		 //# of bytes at left of seed on bseq
		for(j=0; j!=m; j++) {
			_hit.chr=_refid[j];
			_hit.loc=(_refloc[j]<<2)-profile[1][i].a+i;
			if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
				continue;
			//special for mRNA tag alignment
			if((param.tag_type!=-1) && ((_hit.loc<4) || UnequalTag_0(_hit.chr, _hit.loc, ref)))
				continue;			
			d=_hit.loc%12;
			w=CountMismatch(bseq[d], reg[d], ref.bfa[_hit.chr].s+_hit.loc/12);	
//			cout<<"hit:  "<<_hit.loc<<"  "<<w<<endl;

			if(w>param.max_snp_num)
				continue;	
			if(_cur_n_hit[w]>=param.max_num_hits)
				continue;
			if(_tmp_n_hit[w] &&HitExist(_hit, hits[w], hits[w]+_tmp_n_hit[w]))
				continue;
			_hit.z=d;
			hits[w][_cur_n_hit[w]++]=_hit;
		}
	}
	}
	//complementary chain
	if((0==param.chains) ||(2==param.chains)) {
	for(i=0; i!=param.max_snp_num+1; i++) {
		_tmp_n_hit[i]=_cur_n_chit[i];
		if(_cur_n_chit[i])
			sort(chits[i], chits[i]+_cur_n_chit[i], HitComp);
	}
	for(i=0; i!=4; i++) {
		_seed=cseeds[1][i];
//		cout<<"compl: "<<i<<"  seed:  "<<_seed<<"  "<<ref.index[_seed].n1<<endl;
		if((m=ref.index[_seed].n1) ==0) continue;  //no match
//		cout<<i<<"  seed:  "<<_seed<<"  "<<ref.index[_seed].n1<<endl;
		_refid=ref.index[_seed].id1;
		_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
		h=profile[1][i].b1;		 //# of bytes at left of seed on bseq
		for(j=0; j!=m; j++) {
			_hit.chr=_refid[j];
			_hit.loc=(_refloc[j]<<2)-profile[1][i].a+i;
			if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
				continue;
			//special for mRNA tag alignment
			if((param.tag_type!=-1) && UnequalTag_1(_hit.chr, _hit.loc, ref))
				continue;		
			d=_hit.loc%12;			
			w=CountMismatch(cbseq[d], creg[d], ref.bfa[_hit.chr].s+_hit.loc/12);
//			cout<<"hit:  "<<_hit.loc<<"  "<<w<<endl;			

			if(w>param.max_snp_num)
				continue;				
			if(_cur_n_chit[w]>=param.max_num_hits)
				continue;
			if(_tmp_n_hit[w] &&HitExist(_hit, chits[w], chits[w]+_tmp_n_hit[w]))
				continue;	
			_hit.z=d;
			chits[w][_cur_n_chit[w]++]=_hit;
		}
	}
	}
}

char * StrSeed(bit32_t seed, bit32_t size)
{
	char *s = new char[size+1];
	for(int i=size-1; i>=0; i--) {
		s[size-1-i]=param.useful_nt[(seed>>(i*2))&0x3];
	}
	return s;
};

/*n=0: ab; 1: cd; 2: bc; 3: ac; 4: bd; 5: ad*/
//for seed bc ac bd ad
void SingleAlign::SnpAlign_2(RefSeq &ref)
{
	bit32_t i,j,k,m,g,h,w,d;
	for(k=2; k!=6; k++) {
		//direct chain
		if((0==param.chains) ||(1==param.chains)) {
		for(i=0; i<=param.max_snp_num; i++) {
			_tmp_n_hit[i]=_cur_n_hit[i];
			if(_cur_n_hit[i])
				sort(hits[i], hits[i]+_cur_n_hit[i], HitComp);
		}		
		for(i=0; i!=4; i++) {
			_seed=seeds[k][i];
//			cout<<"seed: "<<param.seed_size<<"  "<<k<<"  "<<i<<"  "<<seeds[k][i]<<"  "<<StrSeed(seeds[k][i], param.seed_size)<<endl;
			if(k<=2) {
				if((m=ref.index[_seed].n1) ==0) continue;  //no match
				_refid=ref.index[_seed].id1;
				_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
			}
			else if(k<=4) {
				if((m=ref.index[_seed].n2) ==0) continue;  //no match
				_refid=ref.index[_seed].id2;
				_refloc=ref.index[_seed].loc2;  //list of seeded reference locations				
			}
			else {
				if((m=ref.index[_seed].n3) ==0) continue;  //no match
				_refid=ref.index[_seed].id3;
				_refloc=ref.index[_seed].loc3;  //list of seeded reference locations
			}				
//			cout<<i<<"  seed:  "<<_seed<<"  "<<m<<endl;
			h=profile[k][i].b1;		 //# of bytes at left of seed on bseq
			for(j=0; j!=m; j++) {
				_hit.chr=_refid[j];
				_hit.loc=(_refloc[j]<<2)-profile[k][i].a+i;		
				if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
					continue;		
				//special for mRNA tag alignment
				if((param.tag_type!=-1) && ((_hit.loc<4) || UnequalTag_0(_hit.chr, _hit.loc, ref)))
					continue;				
				d=_hit.loc%12;						
				w=CountMismatch(bseq[d], reg[d], ref.bfa[_hit.chr].s+_hit.loc/12);
//				cout<<"hit:  "<<_hit.loc<<"  "<<w<<endl;

				if(w>param.max_snp_num)
					continue;				
				if(_cur_n_hit[w]>=param.max_num_hits)
					continue;
				if(_tmp_n_hit[w] &&HitExist(_hit, hits[w], hits[w]+_tmp_n_hit[w]))
					continue;					
				_hit.z=d;
				hits[w][_cur_n_hit[w]++]=_hit;
			}
		}
		}
		//complementary chain
		if((0==param.chains) ||(2==param.chains)) {
		for(i=0; i<=param.max_snp_num; i++) {
			_tmp_n_hit[i]=_cur_n_chit[i];
			if(_cur_n_chit[i])
				sort(chits[i], chits[i]+_cur_n_chit[i], HitComp);
		}		
		for(i=0; i!=4; i++) {
			_seed=cseeds[k][i];
			if(k<=2) {
				if((m=ref.index[_seed].n1) ==0) continue;  //no match
				_refid=ref.index[_seed].id1;
				_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
			}
			else if(k<=4) {
				if((m=ref.index[_seed].n2) ==0) continue;  //no match
				_refid=ref.index[_seed].id2;
				_refloc=ref.index[_seed].loc2;  //list of seeded reference locations				
			}
			else {
				if((m=ref.index[_seed].n3) ==0) continue;  //no match
				_refid=ref.index[_seed].id3;
				_refloc=ref.index[_seed].loc3;  //list of seeded reference locations
			}
//			cout<<i<<"  seed:  "<<_seed<<"  "<<m<<endl;
			
			h=profile[k][i].b1;		 //# of bytes at left of seed on bseq
			for(j=0; j!=m; j++) {
				_hit.chr=_refid[j];
				_hit.loc=(_refloc[j]<<2)-profile[k][i].a+i;
				if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
					continue;
				//special for mRNA tag alignment
				if((param.tag_type!=-1) && UnequalTag_1(_hit.chr, _hit.loc, ref))
					continue;		
				d=_hit.loc%12;
				w=CountMismatch(cbseq[d], creg[d], ref.bfa[_hit.chr].s+_hit.loc/12);
//				cout<<"hit:  "<<_hit.loc<<"  "<<w<<endl;

				if(w>param.max_snp_num)
					continue;				
				if(_cur_n_chit[w]>=param.max_num_hits)
					continue;
				if(_tmp_n_hit[w] &&HitExist(_hit, chits[w], chits[w]+_tmp_n_hit[w]))
					continue;		
				_hit.z=d;
				chits[w][_cur_n_chit[w]++]=_hit;
			}
		}
		}
	}	
}

//chain: 0: direct; 1: complementary
int SingleAlign::SnpAlign_range(bool chain, ref_id_t id, ref_loc_t left_end, ref_loc_t right_end, RefSeq &ref)
{
	bit32_t i,j,m,k,g,h,w,z,d;
	ref_id_t *middle_id;
	ref_id_t *end_id;
	ref_loc_t *middle_loc;
	ref_loc_t *end_loc;
	ref_id_t *final_end_id;
	for(i=0; i!=param.max_snp_num+1; i++)
		_cur_n_boundhit[i]=0;
	if(!chain) {   //direct
		for(k=0; k!=6; k++) {
			for(i=0; i<=param.max_snp_num; i++) {
				_tmp_n_hit[i]=_cur_n_boundhit[i];
				if(_cur_n_boundhit[i])
					sort(bound_hits[i], bound_hits[i]+_cur_n_boundhit[i], HitComp);
			}			
			for(i=0; i!=4; i++) {
				_lb=(left_end+profile[k][i].a+i)>>2;
				_rb=(right_end+profile[k][i].a+i)>>2;
				_seed=seeds[k][i];
				if(k<=2) {
					if((m=ref.index[_seed].n1) ==0) continue;  //no match
					_refid=ref.index[_seed].id1;
					_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
				}
				else if(k<=4) {
					if((m=ref.index[_seed].n2) ==0) continue;  //no match
					_refid=ref.index[_seed].id2;
					_refloc=ref.index[_seed].loc2;  //list of seeded reference locations				
				}
				else {
					if((m=ref.index[_seed].n3) ==0) continue;  //no match
					_refid=ref.index[_seed].id3;
					_refloc=ref.index[_seed].loc3;  //list of seeded reference locations
				}
				//	
				end_id=_refid+m-1;
				end_loc=_refloc+m-1;
				final_end_id=_refid+m;
				while(end_id>_refid) {
					middle_id=_refid+(end_id-_refid)/2;
					middle_loc=_refloc+(end_loc-_refloc)/2;
					if(*middle_id<id) {
						_refid=middle_id+1;
						_refloc=middle_loc+1;
					}
					else if(*middle_id>id) {
						end_id=middle_id;
						end_loc=middle_loc;
					}
					else if(*middle_loc<_lb) {
						_refid=middle_id+1;
						_refloc=middle_loc+1;
					}
					else {
						end_id=middle_id;
						end_loc=middle_loc;
					}
				}
				if((*_refid!=id)||(*_refloc<_lb)||(*_refloc>_rb))
					continue;	
				for(; (*_refid==id)&&(_refid<final_end_id)&&(*_refloc<=_rb); _refid++,_refloc++) {	
					_hit.chr=*_refid;
					_hit.loc=(*_refloc<<2)-profile[k][i].a+i;
					if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
						continue;		
					d=_hit.loc%12;
					w=CountMismatch(bseq[d], reg[d], ref.bfa[_hit.chr].s+_hit.loc/12);
					if(w >param.max_snp_num)
						continue;
					if(HitExist(_hit, bound_hits[w], bound_hits[w]+_tmp_n_hit[w]))
						continue;			
					_hit.z=d;
					bound_hits[w][_cur_n_boundhit[w]++]=_hit;
				}	
			}
		}
	}
	else {  //complementary
		for(k=0; k!=6; k++) {
			for(i=0; i<=param.max_snp_num; i++) {
				_tmp_n_hit[i]=_cur_n_boundhit[i];
				if(_cur_n_boundhit[i])
					sort(bound_hits[i], bound_hits[i]+_cur_n_boundhit[i], HitComp);
			}			
			for(i=0; i!=4; i++) {
				_lb=(left_end+profile[k][i].a+i)>>2;
				_rb=(right_end+profile[k][i].a+i)>>2;
				//so the other end should be on complementary chain
				_seed=cseeds[k][i];
				if(k<=2) {
					if((m=ref.index[_seed].n1) ==0) continue;  //no match
					_refid=ref.index[_seed].id1;
					_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
				}
				else if(k<=4) {
					if((m=ref.index[_seed].n2) ==0) continue;  //no match
					_refid=ref.index[_seed].id2;
					_refloc=ref.index[_seed].loc2;  //list of seeded reference locations				
				}
				else {
					if((m=ref.index[_seed].n3) ==0) continue;  //no match
					_refid=ref.index[_seed].id3;
					_refloc=ref.index[_seed].loc3;  //list of seeded reference locations
				}
				//
				end_id=_refid+m-1;
				end_loc=_refloc+m-1;
				final_end_id=_refid+m;
				while(end_id>_refid) {
					middle_id=_refid+(end_id-_refid)/2;
					middle_loc=_refloc+(end_loc-_refloc)/2;
					if(*middle_id<id) {
						_refid=middle_id+1;
						_refloc=middle_loc+1;
					}
					else if(*middle_id>id) {
						end_id=middle_id;
						end_loc=middle_loc;
					}
					else if(*middle_loc<_lb) {
						_refid=middle_id+1;
						_refloc=middle_loc+1;
					}
					else {
						end_id=middle_id;
						end_loc=middle_loc;
					}			
				}
				if((*_refid!=id)||(*_refloc<_lb)||(*_refloc>_rb))
					continue;	
				for(; (*_refid==id)&&(_refid<final_end_id)&&(*_refloc<=_rb); _refid++,_refloc++) {		
					_hit.chr=*_refid;
					_hit.loc=(*_refloc<<2)-profile[k][i].a+i;
					if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()) //overflow the end of refseq
						continue;			
					d=_hit.loc%12;
					w=CountMismatch(cbseq[d], creg[d], ref.bfa[_hit.chr].s+_hit.loc/12);
					if(w >param.max_snp_num)
						continue;
					if(HitExist(_hit, bound_hits[w], bound_hits[w]+_tmp_n_hit[w]))
						continue;								
					_hit.z=d;
					bound_hits[w][_cur_n_boundhit[w]++]=_hit;
				}
			}
		}
	}	
	for(i=0; i<=param.max_snp_num; i++) {
		if(_cur_n_boundhit[i])
			return i;
	}
	return -1;
}

//delete g-bps starting with pos
void SingleAlign::DelBp(bit24_t *ori, bit24_t *ori_reg, bit32_t pos, bit32_t g)
{
	int i;
	int n=pos/12;
	int k=pos%12;
	//pre elements
	for(i=0; i<n; i++) {
		gbseq[i]=ori[i];
		greg[i]=ori_reg[i];
	}
	if(i<FIXELEMENT-1) {
  	//the element
  	if(k+g>=12) {
  		gbseq[i].a=(ori[i].a&_leftbits[k])|((ori[i+1].a>>(24-2*g))&_rightbits[12-k]);
  		greg[i].a=(ori_reg[i].a&_leftbits[k])|((ori_reg[i+1].a>>(24-2*g))&_rightbits[12-k]);
  	}
  	else {
  		gbseq[i].a=(ori[i].a&_leftbits[k])|(((ori[i].a<<(2*g))|(ori[i+1].a>>(24-2*g)))&_rightbits[12-k]);
  		greg[i].a=(ori_reg[i].a&_leftbits[k])|(((ori_reg[i].a<<(2*g))|(ori_reg[i+1].a>>(24-2*g)))&_rightbits[12-k]);
  	}
  	i++;
  	//following elements
  	for(; i<FIXELEMENT-1; i++) {
  			gbseq[i].a=(ori[i].a<<(2*g))|(ori[i+1].a>>(24-2*g));
  			greg[i].a=(ori_reg[i].a<<(2*g))|(ori_reg[i+1].a>>(24-2*g));
  	}
  	//last element
  	gbseq[i].a=ori[i].a<<(2*g);
  	greg[i].a=ori_reg[i].a<<(2*g);
  }
  else {
  	//the element
  	if(k+g>=12) {
  		gbseq[i].a=ori[i].a&_leftbits[k];
  		greg[i].a=ori_reg[i].a&_leftbits[k];
  	}
  	else {
  		gbseq[i].a=(ori[i].a&_leftbits[k])|((ori[i].a<<(2*g))&_rightbits[12-k]);
  		greg[i].a=(ori_reg[i].a&_leftbits[k])|((ori_reg[i].a<<(2*g))&_rightbits[12-k]);
  	}  	
  }
}

//insert g-bps in front of pos
void SingleAlign::InsBp(bit24_t *ori, bit24_t *ori_reg, bit32_t pos, bit32_t g)
{
	int i;
	int n=pos/12;
	int k=pos%12;
//	cout<<pos<<"  "<<g<<"  n: "<<n<<"  k: "<<k<<endl;
	//pre elements
	for(i=0; i<n; i++) {
		gbseq[i]=ori[i];
		greg[i]=ori_reg[i];
	}
	//the element
	if(k+g>=12) {
		gbseq[i].a=ori[i].a&_leftbits[k];
		greg[i].a=ori_reg[i].a&_leftbits[k];
		i++;
		gbseq[i].a=((ori[i-1].a&_rightbits[12-k])<<(24-2*g))|(ori[i].a>>(2*g));
		greg[i].a=((ori_reg[i-1].a&_rightbits[12-k])<<(24-2*g))|(ori_reg[i].a>>(2*g));
	}
	else {
		gbseq[i].a=(ori[i].a&_leftbits[k])|((ori[i].a>>(2*g))&_rightbits[12-(k+g)]);
		greg[i].a=(ori_reg[i].a&_leftbits[k])|((ori_reg[i].a>>(2*g))&_rightbits[12-(k+g)]);
	}
	i++;
	//following elements
	for(; i<FIXELEMENT; i++) {
		gbseq[i].a=((ori[i-1].a&_rightbits[g])<<(24-2*g))|(ori[i].a>>(2*g));
		greg[i].a=((ori_reg[i-1].a&_rightbits[g])<<(24-2*g))|(ori_reg[i].a>>(2*g));
	}
	//last element
	gbseq[i].a=(ori[i-1].a&_rightbits[g])<<(24-2*g);
	greg[i].a=(ori_reg[i-1].a&_rightbits[g])<<(24-2*g);
/*
	for(i=0; i<=FIXELEMENT; i++) {
		cout<<greg[i].a<<"  ";
	}
	cout<<endl;
*/
}

int SingleAlign::GapAlign(RefSeq &ref)
{
	bit32_t i,j,m,k,g,h,w,z,d;
	/* break as soon as get smaller indels;
	   break as soon as get left most gapped align */
	/*
			  ==============================
		  0 |_____                 0_____| 
		  1  |_____               1_____|
		  2   |______            2_____|
		  3    |______          3_____|
	*/
	
	//make seeds and get candidates
	for(i=0; i<4; i++) {
		_gc1[0][i]= (i+param.seed_size<=12)? (bseq[0][0].a>>2*(12-i-param.seed_size))&param.seed_bits : (bseq[0][0].a<<2*(param.seed_size-(12-i))|bseq[0][1].a>>2*(24-i-param.seed_size))&param.seed_bits;
		_gc2[0][i]= (i+param.seed_size<=12)? (cbseq[0][0].a>>2*(12-i-param.seed_size))&param.seed_bits : (cbseq[0][0].a<<2*(param.seed_size-(12-i))|cbseq[0][1].a>>2*(24-i-param.seed_size))&param.seed_bits;
		j=(_pread->seq.size()-param.seed_size-i)%12;
		k=(_pread->seq.size()-param.seed_size-i)/12;
		_gc1[1][i]= (j+param.seed_size<=12)? (bseq[0][k].a>>2*(12-j-param.seed_size))&param.seed_bits : (bseq[0][k].a<<2*(param.seed_size-(12-j))|bseq[0][k+1].a>>2*(24-j-param.seed_size))&param.seed_bits;
		_gc2[1][i]= (j+param.seed_size<=12)? (cbseq[0][k].a>>2*(12-j-param.seed_size))&param.seed_bits : (cbseq[0][k].a<<2*(param.seed_size-(12-j))|cbseq[0][k+1].a>>2*(24-j-param.seed_size))&param.seed_bits;
		
		_num1[0][i]=ref.index[_gc1[0][i]].n1;
		_num1[1][i]=ref.index[_gc1[1][i]].n1;
		_num2[0][i]=ref.index[_gc2[0][i]].n1;
		_num2[1][i]=ref.index[_gc2[1][i]].n1;
		_id1[0][i]=ref.index[_gc1[0][i]].id1;
		_loc1[0][i]=ref.index[_gc1[0][i]].loc1;
		_id1[1][i]=ref.index[_gc1[1][i]].id1;
		_loc1[1][i]=ref.index[_gc1[1][i]].loc1;	
		_id2[0][i]=ref.index[_gc2[0][i]].id1;
		_loc2[0][i]=ref.index[_gc2[0][i]].loc1;
		_id2[1][i]=ref.index[_gc2[1][i]].id1;
		_loc2[1][i]=ref.index[_gc2[1][i]].loc1;
//		cout<<i<<"  "<<_num1[0][i]<<"  "<<_num1[1][i]<<"  "<<_num2[0][i]<<"  "<<_num2[1][i]<<endl;
//		for(m=0; m<_num1[0][i]; m++)
//			cout<<_loc1[0][i][m]*4<<endl;
	}
	//try alignment by insert or delete gaps
	for(g=1; g<=param.max_gap_size; g++) {
		_gap_size=g;
		for(i=0; i<4; i++) {
			//direct chain
			if((0==param.chains) ||(1==param.chains)) {
			//left part of seeds
			for(m=0, _tid=_id1[0][i], _tp=_loc1[0][i]; m<_num1[0][i]; m++, _tid++, _tp++) {
				_hit.chr=*_tid;
				_hit.loc=(*_tp<<2)-i;
				if((param.tag_type!=-1) && (UnequalTag_0(_hit.chr, _hit.loc, ref)))
					continue;						
				if(ref.title[_hit.chr].size<_hit.loc+_pread->seq.size()+g)
					continue;
				k=_hit.loc%12;
				h=_hit.loc/12;				
				//insertion on read, so try to delete the bps and align
				if(_cur_n_gaphit >= param.max_num_hits)
					break;				
				for(j=_pread->seq.size()-param.gap_edge-g; j>=i+param.seed_size; j--) {
					DelBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif	
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a)||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif		
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif		
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif		
					_hit.z=100+j;
					gaphits[_cur_n_gaphit++]=_hit;
					break;
				}
				//deletion on read, so try to insert bps in read and align
				if(_cur_n_gaphit >= param.max_num_hits)
					break;					
				for(j=_pread->seq.size()-param.gap_edge; j>=i+param.seed_size; j--) {
//					cout<<"gap:  "<<g<<"  coord:  "<<_hit.loc<<"  offset: "<<j<<"  "<<k<<endl;
					InsBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
/*
					for(int mmm=0; mmm<=6; mmm++)
						cout<<ref.bfa[_hit.chr].s[h+mmm].a<<"  "<<gbseq[mmm].a<<"  "<<greg[mmm].a<<"  "
						<<((ref.bfa[_hit.chr].s[h+mmm].a^gbseq[mmm].a)&greg[mmm].a)<<endl;
*/
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a)  ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					gaphits[_cur_n_gaphit++]=_hit;
					break;
				}
			}
			sort(gaphits, gaphits+_cur_n_gaphit, HitComp);
			_tmp_n_gaphit=_cur_n_gaphit;
			//right part of seeds
			for(m=0, _tid=_id1[1][i], _tp=_loc1[1][i]; m<_num1[1][i]; m++, _tid++, _tp++) {		
				//insertion on read, so try to delete the bps and align
				if(_cur_n_gaphit >= param.max_num_hits)
					break;
				_hit.chr=*_tid;
				_hit.loc=(*_tp<<2)-(_pread->seq.size()-param.seed_size-i-g);
				if(ref.title[_hit.chr].size<_hit.loc+_pread->seq.size()+g)
					continue;
				if(HitExist(_hit, gaphits, gaphits+_tmp_n_gaphit))
					continue;
				k=_hit.loc%12;
				h=_hit.loc/12;
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					DelBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif	
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif	
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif	
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif	
					_hit.z=100+j;
					gaphits[_cur_n_gaphit++]=_hit;
					break;
				}
				//deletion on read
				if(_cur_n_gaphit >= param.max_num_hits)
					break;	
				_hit.chr=*_tid;
				_hit.loc=(*_tp<<2)-(_pread->seq.size()-param.seed_size-i+g);
				if(ref.title[_hit.chr].size<_hit.loc+_pread->seq.size()+g)
					continue;					
				if(HitExist(_hit, gaphits, gaphits+_tmp_n_gaphit))
					continue;					
				k=_hit.loc%12;
				h=_hit.loc/12;									
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					InsBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a) ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					gaphits[_cur_n_gaphit++]=_hit;
					break;
				}
			}		
			}	
			//reverse chain
			if((0==param.chains) ||(-1==param.chains)) {
			//left part of seeds
			for(m=0, _tid=_id2[0][i], _tp=_loc2[0][i]; m<_num2[0][i]; m++, _tid++, _tp++) {
				_hit.chr=*_tid;
				_hit.loc=(*_tp<<2)-i;
				if((param.tag_type!=-1) && (UnequalTag_1(_hit.chr, _hit.loc, ref)))
					continue;						
				if(ref.title[_hit.chr].size<_hit.loc+_pread->seq.size()+g)
					continue;
				k=_hit.loc%12;
				h=_hit.loc/12;				
				//insertion on read, so try to delete the bps and align
				if(_cur_n_cgaphit >= param.max_num_hits)
					break;				
				for(j=_pread->seq.size()-param.gap_edge-g; j>=i+param.seed_size; j--) {
					DelBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif	
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif	
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif	
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif	
					_hit.z=100+j;
					cgaphits[_cur_n_cgaphit++]=_hit;
					break;
				}
				//deletion on read, so try to insert bps in read and align
				if(_cur_n_cgaphit >= param.max_num_hits)
					break;					
				for(j=_pread->seq.size()-param.gap_edge; j>=i+param.seed_size; j--) {
					InsBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a) ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					cgaphits[_cur_n_cgaphit++]=_hit;
					break;
				}
			}
			sort(cgaphits, cgaphits+_cur_n_cgaphit, HitComp);
			_tmp_n_cgaphit=_cur_n_cgaphit;			
			//right part of seeds
			for(m=0, _tid=_id2[1][i], _tp=_loc2[1][i]; m<_num2[1][i]; m++, _tid++, _tp++) {		
				//insertion on read, so try to delete the bps and align
				if(_cur_n_cgaphit >= param.max_num_hits)
					break;
				_hit.chr=*_tid;
				_hit.loc=(*_tp<<2)-(_pread->seq.size()-param.seed_size-i-g);
				if(ref.title[_hit.chr].size<_hit.loc+_pread->seq.size()+g)
					continue;
				if(HitExist(_hit, cgaphits, cgaphits+_tmp_n_cgaphit))
					continue;					
				k=_hit.loc%12;
				h=_hit.loc/12;
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					DelBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif	
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif	
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) )
						continue;
#endif	
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif	
					_hit.z=100+j;
					cgaphits[_cur_n_cgaphit++]=_hit;
					break;
				}
				//deletion on read
				if(_cur_n_cgaphit >= param.max_num_hits)
					break;	
				_hit.chr=*_tid;
				_hit.loc=(*_tp<<2)-(_pread->seq.size()-param.seed_size-i+g);
				if(ref.title[_hit.chr].size<_hit.loc+_pread->seq.size()+g)
					continue;
				if(HitExist(_hit, cgaphits, cgaphits+_tmp_n_cgaphit))
					continue;						
				k=_hit.loc%12;
				h=_hit.loc/12;									
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					InsBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a) ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					cgaphits[_cur_n_cgaphit++]=_hit;
//					cout<<"ok here2"<<endl;
					break;
				}
			}
			}	
		}
		if(_cur_n_gaphit ||_cur_n_cgaphit)
			return _cur_n_gaphit+_cur_n_cgaphit;
	}
	return 0;
}

//return gap size
//chain: 0: direct; 1: complementary
int SingleAlign::GapAlign_range(bool chain, ref_id_t id, ref_loc_t left_end, ref_loc_t right_end, RefSeq &ref)
{
	_cur_n_boundgaphit=0;
	bit32_t i,j,m,k,g,h,w,z,d,x;
	ref_id_t *middle_id;
	ref_id_t *end_id;
	ref_id_t *final_end_id;
	ref_loc_t *middle_loc;
	ref_loc_t *end_loc;
	/* break as soon as get smaller indels;
	   break as soon as get left most gapped align */
	/*
			  ==============================
		  0 |_____                 0_____| 
		  1  |_____               1_____|
		  2   |______            2_____|
		  3    |______          3_____|
	*/
	if(!chain) {  //direct chain
		//make seeds and get candidates
		for(i=0; i<4; i++) {
			_gc1[0][i]= (i<=12-param.seed_size)? (bseq[0][0].a>>2*(12-i-param.seed_size))&param.seed_bits : (bseq[0][0].a<<2*(param.seed_size-(12-i))|bseq[0][1].a>>2*(24-i-param.seed_size))&param.seed_bits;
			j=(_pread->seq.size()-param.seed_size-i)%12;
			k=(_pread->seq.size()-param.seed_size-i)/12;
			_gc1[1][i]= (12-j>=param.seed_size)? (bseq[0][k].a>>2*(12-j-param.seed_size))&param.seed_bits : (bseq[0][k].a<<2*(param.seed_size-(12-j))|bseq[0][k+1].a>>2*(24-j-param.seed_size))&param.seed_bits;
		
			if(ref.index[_gc1[0][i]].n1) {
				_lb=(left_end+i-param.max_gap_size)>>2;
				_rb=(right_end+i+param.max_gap_size)>>2;
				
  			_refid=ref.index[_gc1[0][i]].id1;
  			_refloc=ref.index[_gc1[0][i]].loc1;
  			end_id=_refid+ref.index[_gc1[0][i]].n1-1;
  			end_loc=_refloc+ref.index[_gc1[0][i]].n1-1;
  			final_end_id=_refid+ref.index[_gc1[0][i]].n1;
  			while(end_id>_refid) {
  				middle_id=_refid+(end_id-_refid)/2;
  				middle_loc=_refloc+(end_loc-_refloc)/2;
  				if(*middle_id<id) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else if(*middle_id>id) {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}
  				else if(*middle_loc<_lb) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}			
  			}
  			if((*_refid !=id) || (*_refloc<_lb) || (*_refloc>_rb))
  				_num1[0][i]=0;
  			else {
  				_num1[0][i]=final_end_id-_refid;
  				_id1[0][i]=_refid;
  				_loc1[0][i]=_refloc;
  			}
  		}
  		else
  			_num1[0][i]=0;
			
			if(ref.index[_gc1[1][i]].n1) {
				_lb=(left_end+_pread->seq.size()-param.seed_size-i-param.max_gap_size)>>2;
				_rb=(right_end+_pread->seq.size()-param.seed_size-i+param.max_gap_size)>>2;
				
  			_refid=ref.index[_gc1[1][i]].id1;
  			_refloc=ref.index[_gc1[1][i]].loc1;
  			end_id=_refid+ref.index[_gc1[1][i]].n1-1;
  			end_loc=_refloc+ref.index[_gc1[1][i]].n1-1;
  			final_end_id=_refid+ref.index[_gc1[1][i]].n1;
  			while(end_id>_refid) {
  				middle_id=_refid+(end_id-_refid)/2;
  				middle_loc=_refloc+(end_loc-_refloc)/2;
  				if(*middle_id<id) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else if(*middle_id>id) {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}
  				else if(*middle_loc<_lb) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}			
  			}
  			if((*_refid !=id) || (*_refloc<_lb) || (*_refloc>_rb))
  				_num1[1][i]=0;
  			else {
  				_num1[1][i]=final_end_id-_refid;
  				_id1[1][i]=_refid;
  				_loc1[1][i]=_refloc;
  			}
  		}
  		else
  			_num1[1][i]=0;
		}	
		//try alignment by insert or delete gaps
		for(g=1; g<=param.max_gap_size; g++) {
			_gap_size=g;
			for(i=0; i<4; i++) {
			//left part of seeds
			_refid=_id1[0][i];
			_refloc=_loc1[0][i];
			for(x=0; (x<_num1[0][i])&&(*_refid==id); _refid++,_refloc++,x++) {
				_hit.chr=*_refid;
				_hit.loc=(*_refloc<<2)-i;
				if(_hit.loc>right_end)
					break;
				if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()-g)
					continue;
				k=_hit.loc%12;
				h=_hit.loc/12;				
				//insertion on read, so try to delete the bps and align
				if(_cur_n_gaphit >= param.max_num_hits)
					break;				
				for(j=_pread->seq.size()-param.gap_edge-g; j>=i+param.seed_size; j--) {
					DelBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif	
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a)||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif		
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif		
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif		
					_hit.z=100+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
				//deletion on read, so try to insert bps in read and align
				if(_cur_n_boundgaphit >= param.max_num_hits)
					break;					
				for(j=_pread->seq.size()-param.gap_edge; j>=i+param.seed_size; j--) {
					InsBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a) ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
			}
			sort(bound_gaphits, bound_gaphits+_cur_n_boundgaphit, HitComp);
			_tmp_n_gaphit=_cur_n_boundgaphit;			
			//right part of seeds
			_refid=_id1[1][i];
			_refloc=_loc1[1][i];
			for(x=0; (x<_num1[1][i])&&(*_refid==id); _refid++,_refloc++,x++) {
				//insertion on read, so try to delete the bps and align
				if(_cur_n_boundgaphit >= param.max_num_hits)
					break;
				_hit.chr=*_refid;
				_hit.loc=(*_refloc<<2)-(_pread->seq.size()-param.seed_size-i-g);
				if(_hit.loc>right_end)
					break;
				if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()-g)
					continue;
				if(HitExist(_hit, bound_gaphits, bound_gaphits+_tmp_n_gaphit))
					continue;
				k=_hit.loc%12;
				h=_hit.loc/12;
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					DelBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif	
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif	
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif	
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif	
					_hit.z=100+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
				//deletion on read
				if(_cur_n_gaphit >= param.max_num_hits)
					break;	
				_hit.chr=*_refid;
				_hit.loc=(*_refloc<<2)-(_pread->seq.size()-param.seed_size-i+g);
				if(_hit.loc>right_end)
					break;				
				if(ref.title[_hit.chr].size<_hit.loc+_pread->seq.size()+g)
					continue;
				if(HitExist(_hit, bound_gaphits, bound_gaphits+_tmp_n_gaphit))
					continue;					
				k=_hit.loc%12;
				h=_hit.loc/12;									
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					InsBp(bseq[k], reg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a) ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
			}		
			}
			if(_cur_n_boundgaphit)
				return g;			
		}
	}
	
	else {  //complementary
		//make seeds and get candidates
		for(i=0; i<4; i++) {
			_gc2[0][i]= (i<=12-param.seed_size)? (cbseq[0][0].a>>2*(12-i-param.seed_size))&param.seed_bits : (cbseq[0][0].a<<2*(param.seed_size-(12-i))|cbseq[0][1].a>>2*(24-i-param.seed_size))&param.seed_bits;
			j=(_pread->seq.size()-param.seed_size-i)%12;
			k=(_pread->seq.size()-param.seed_size-i)/12;
			_gc2[1][i]= (12-j>=param.seed_size)? (cbseq[0][k].a>>2*(12-j-param.seed_size))&param.seed_bits : (cbseq[0][k].a<<2*(param.seed_size-(12-j))|cbseq[0][k+1].a>>2*(24-j-param.seed_size))&param.seed_bits;
		
			if(ref.index[_gc2[0][i]].n1) {
				_lb=(left_end+i-param.max_gap_size)>>2;
				_rb=(right_end+i+param.max_gap_size)>>2;
								
  			_refid=ref.index[_gc2[0][i]].id1;
  			_refloc=ref.index[_gc2[0][i]].loc1;
  			end_id=_refid+ref.index[_gc2[0][i]].n1-1;
  			end_loc=_refloc+ref.index[_gc2[0][i]].n1-1;
  			final_end_id=_refid+ref.index[_gc2[0][i]].n1;
  			while(end_id>_refid) {
  				middle_id=_refid+(end_id-_refid)/2;
  				middle_loc=_refloc+(end_loc-_refloc)/2;
  				if(*middle_id<id) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else if(*middle_id>id) {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}
  				else if(*middle_loc<_lb) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}			
  			}
  			if((*_refid !=id) || (*_refloc<_lb) || (*_refloc>_rb))
  				_num2[0][i]=0;
  			else {
  				_num2[0][i]=final_end_id-_refid;
  				_id2[0][i]=_refid;
  				_loc2[0][i]=_refloc;
  			}
  		}
  		else
  			_num2[0][i]=0;
			
			if(ref.index[_gc2[1][i]].n1) {
				_lb=(left_end+_pread->seq.size()-param.seed_size-i-param.max_gap_size)>>2;
				_rb=(right_end+_pread->seq.size()-param.seed_size-i+param.max_gap_size)>>2;
								
  			_refid=ref.index[_gc2[1][i]].id1;
  			_refloc=ref.index[_gc2[1][i]].loc1;
  			end_id=_refid+ref.index[_gc2[1][i]].n1-1;
  			end_loc=_refloc+ref.index[_gc2[1][i]].n1-1;
  			final_end_id=_refid+ref.index[_gc2[1][i]].n1;
  			while(end_id>_refid) {
  				middle_id=_refid+(end_id-_refid)/2;
  				middle_loc=_refloc+(end_loc-_refloc)/2;
  				if(*middle_id<id) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else if(*middle_id>id) {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}
  				else if(*middle_loc<_lb) {
  					_refid=middle_id+1;
  					_refloc=middle_loc+1;
  				}
  				else {
  					end_id=middle_id;
  					end_loc=middle_loc;
  				}			
  			}
  			if((*_refid !=id) || (*_refloc<_lb) || (*_refloc>_rb))
  				_num2[1][i]=0;
  			else {
  				_num2[1][i]=final_end_id-_refid;
  				_id2[1][i]=_refid;
  				_loc2[1][i]=_refloc;
  			}
  		}
  		else
  			_num2[1][i]=0;
		}	
		//try alignment by insert or delete gaps
		for(g=1; g<=param.max_gap_size; g++) {
			_gap_size=g;
			for(i=0; i<4; i++) {
			//left part of seeds
			_refid=_id2[0][i];
			_refloc=_loc2[0][i];
			for(x=0; (x<_num2[0][i])&&(*_refid==id); _refid++,_refloc++,x++) {
				_hit.chr=*_refid;
				_hit.loc=(*_refloc<<2)-i;
				if(_hit.loc>right_end)
					break;				
				if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()-g)
					continue;
				k=_hit.loc%12;
				h=_hit.loc/12;				
				//insertion on read, so try to delete the bps and align
				if(_cur_n_boundgaphit >= param.max_num_hits)
					break;				
				for(j=_pread->seq.size()-param.gap_edge-g; j>=i+param.seed_size; j--) {
					DelBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif	
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif	
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif	
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif	
					_hit.z=100+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
				//deletion on read, so try to insert bps in read and align
				if(_cur_n_boundgaphit >= param.max_num_hits)
					break;					
				for(j=_pread->seq.size()-param.gap_edge; j>=i+param.seed_size; j--) {
					InsBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a) ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
			}
			sort(bound_gaphits, bound_gaphits+_cur_n_boundgaphit, HitComp);
			_tmp_n_gaphit=_cur_n_boundgaphit;				
			//right part of seeds
			_refid=_id2[1][i];
			_refloc=_loc2[1][i];
			for(x=0; (x<_num2[1][i])&&(*_refid==id); _refid++,_refloc++,x++) {	
				//insertion on read, so try to delete the bps and align
				if(_cur_n_boundgaphit >= param.max_num_hits)
					break;
				_hit.chr=*_refid;
				_hit.loc=(*_refloc<<2)-(_pread->seq.size()-param.seed_size-i-g);
				if(_hit.loc>right_end)
					break;				
				if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()-g)
					continue;
				if(HitExist(_hit, bound_gaphits, bound_gaphits+_tmp_n_gaphit))
					continue;
				k=_hit.loc%12;
				h=_hit.loc/12;
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					DelBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a))
						continue;
#endif	
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif	
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif	
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a) ||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif	
					_hit.z=100+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
				//deletion on read
				if(_cur_n_boundgaphit >= param.max_num_hits)
					break;	
				_hit.chr=*_refid;
				_hit.loc=(*_refloc<<2)-(_pread->seq.size()-param.seed_size-i+g);
				if(_hit.loc>right_end)
					break;				
				if(_hit.loc>ref.title[_hit.chr].size-_pread->seq.size()-g)
					continue;
				if(HitExist(_hit, bound_gaphits, bound_gaphits+_tmp_n_gaphit))
					continue;										
				k=_hit.loc%12;
				h=_hit.loc/12;									
				for(j=param.seed_size+i-1; j>=param.gap_edge; j--) {
					InsBp(cbseq[k], creg[k], j+k, g);
#ifdef READ_36
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a))
						continue;
#endif
#ifdef READ_48
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a))
						continue;
#endif
#ifdef READ_60
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a))
						continue;
#endif
#ifdef READ_72
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a))
						continue;
#endif
#ifdef READ_84
					if(((ref.bfa[_hit.chr].s[h].a^gbseq[0].a)&greg[0].a) ||((ref.bfa[_hit.chr].s[h+1].a^gbseq[1].a)&greg[1].a)
						||((ref.bfa[_hit.chr].s[h+2].a^gbseq[2].a)&greg[2].a) ||((ref.bfa[_hit.chr].s[h+3].a^gbseq[3].a)&greg[3].a)
						||((ref.bfa[_hit.chr].s[h+4].a^gbseq[4].a)&greg[4].a) ||((ref.bfa[_hit.chr].s[h+5].a^gbseq[5].a)&greg[5].a)
						||((ref.bfa[_hit.chr].s[h+6].a^gbseq[6].a)&greg[6].a) ||((ref.bfa[_hit.chr].s[h+7].a^gbseq[7].a)&greg[7].a) ||((ref.bfa[_hit.chr].s[h+8].a^gbseq[8].a)&greg[8].a))
						continue;
#endif
					_hit.z=200+j;
					bound_gaphits[_cur_n_boundgaphit++]=_hit;
					break;
				}
			}
			}	
			if(_cur_n_boundgaphit)
				return g;
		}
	}
	return -1;
}

void SingleAlign::ClearHits()
{
	for(int i=0; i<=param.max_snp_num; i++)
		_cur_n_hit[i]=_cur_n_chit[i]=0;
	_cur_n_gaphit=_cur_n_cgaphit=0;	
}
int SingleAlign::RunAlign(RefSeq &ref)
{
	ClearHits();
	ConvertBinaySeq();
	//ab, snp/exact alignment
	GenerateSeeds_1(0);
	SnpAlign_0(ref);
	if(_cur_n_hit[0]||_cur_n_chit[0]) 
		return 1; //end if find exactly matched hits
	
	if(param.max_snp_num>0) {
		GenerateSeeds_1(1);
		SnpAlign_1(ref);
		if(_cur_n_hit[1]||_cur_n_chit[1]||_cur_n_hit[0]||_cur_n_chit[0])
			return 1; //end if find 1 mismatch hits
	}
	
	if(param.max_snp_num>1) {
		//complete snp alignment
		GenerateSeeds_1(2);
		GenerateSeeds_2(3);
		GenerateSeeds_2(5);
		GenerateSeeds_3(4);
		SnpAlign_2(ref);
		for(int i=0; i<=param.max_snp_num; i++) {
			if(_cur_n_hit[i]||_cur_n_chit[i]) 
				return 1;
		}
	}
/*	for(int k=0; k<6; k++)
		for(int i=0; i<4; i++)
			cout<<"seed: "<<param.seed_size<<"  "<<k<<"  "<<i<<"  "<<seeds[k][i]<<"  "<<StrSeed(seeds[k][i], param.seed_size)<<endl;
*/
	//gapped alignment
	if(param.max_gap_size)
		if(GapAlign(ref))
			return 1;
	return 0;
}

int SingleAlign::FilterReads()
{
	if((_pread->seq.size()<param.min_read_size)||(CountNs()>param.max_ns)) {  //filter if too many 'N's
		return 1;
	}
	return 0;
}

int SingleAlign::CountStringMismatch(int offset, string &s1, string s2)
{
	int num_mismatch=0;
	for(int i=0; (i+offset<s1.size())&&(i<s2.size()); i++) {
		if(s1[i+offset]^s2[i]) {
			num_mismatch++;
			if(num_mismatch>param.admis)
				return 0;
		}
	}
	return 1;
}

void SingleAlign::Do_Batch(RefSeq &ref)
{
	_str_align.clear();
	bit32_t tt;
	//for mRNA tag alignment
	if(param.tag_type!=-1) {
		int tag_len=param.tag_seq.size();
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			_pread->seq.erase(0, tag_len);
			_pread->qual.erase(0, tag_len);
			_pread->seq.resize(param.tag_remain);
			_pread->qual.resize(param.tag_remain);
			if((_pread->seq.size()<param.min_read_size)||(CountNs()>param.max_ns)) {
				continue;
			}
			if(RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
			}		
		}
		return;
	}
	//for miRNA alignment
	if(param.adapter.size()) {
		string ori_seq;
		string ori_qual;
		int i;
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			ori_seq=_pread->seq;
			ori_qual=_pread->qual;
			for(i=param.mirna_max; i>=param.mirna_min; i--) {
				if(CountStringMismatch(i, ori_seq, param.adapter)) {
					_pread->seq=ori_seq;
					_pread->qual=ori_qual;
					_pread->seq.erase(i);
					_pread->qual.erase(i);
					if((CountNs()<=param.max_ns)&&RunAlign(ref)) {
						StringAlign(ref, _str_align);
						n_aligned++;
						break;
					}
				}
			}
		}
		return;
	}
	//alignment for normal sequences
	if(0==param.trim_lowQ) {
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			if((_pread->seq.size()<param.min_read_size)||(CountNs()>param.max_ns)) {
				continue;
			}
			if(RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
			}
		}
	}
	else if(10>=param.trim_lowQ) {
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			_pread->seq.erase(_pread->seq.size()-param.trim_lowQ, param.trim_lowQ);
			_pread->qual.erase(_pread->qual.size()-param.trim_lowQ, param.trim_lowQ);
			if((_pread->seq.size()<param.min_read_size)||(CountNs()>param.max_ns)) {
				continue;
			}
			if(RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
			}
		}
	}
	else if(20>=param.trim_lowQ) {
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			_pread->seq.erase(_pread->seq.size()-(param.trim_lowQ-10), param.trim_lowQ-10);
			_pread->qual.erase(_pread->qual.size()-(param.trim_lowQ-10), param.trim_lowQ-10);
			_pread->seq.erase(0,1);
			_pread->qual.erase(0,1);
			if((_pread->seq.size()<param.min_read_size)||(CountNs()>param.max_ns)) {
				continue;
			}
			if(RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
			}
		}		
	}
	else if(30>=param.trim_lowQ) {
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			if((_pread->seq.size()>=param.min_read_size)&&(CountNs()<=param.max_ns)&&RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
			}			
			else {
				_pread->seq.erase(_pread->seq.size()-(param.trim_lowQ-20), param.trim_lowQ-20);
				_pread->qual.erase(_pread->qual.size()-(param.trim_lowQ-20), param.trim_lowQ-20);		
				if((_pread->seq.size()>=param.min_read_size)||(CountNs()<=param.max_ns)) {
					continue;
				}
				if(RunAlign(ref)) {
					StringAlign(ref, _str_align);	
					n_aligned++;
				}			
			}
		}				
	}
	else if(40>=param.trim_lowQ) {
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			if((_pread->seq.size()>=param.min_read_size)&&(CountNs()<=param.max_ns)&&RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
			}			
			else {
				_pread->seq.erase(_pread->seq.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
				_pread->qual.erase(_pread->qual.size()-(param.trim_lowQ-30), param.trim_lowQ-30);
				_pread->seq.erase(0,1);
				_pread->qual.erase(0,1);		
				if((_pread->seq.size()<param.min_read_size)||(CountNs()>param.max_ns)) {
					continue;
				}
				if(RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned++;
				}
			}
		}			
	}
	else if(50>=param.trim_lowQ) {
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			if((_pread->seq.size()>=param.min_read_size)&&(CountNs()<=param.max_ns)&&RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
				continue;
			}	
			while(1) {
				_pread->seq.erase(_pread->seq.size()-(param.trim_lowQ-40), param.trim_lowQ-40);
				_pread->qual.erase(_pread->qual.size()-(param.trim_lowQ-40), param.trim_lowQ-40);						
				if(_pread->seq.size()<param.min_read_size)
					break;
				if((CountNs()<=param.max_ns)&&RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned++;
					break;
				}
			}
		}		
	}
	else if(60>=param.trim_lowQ) {
		for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
			if((_pread->seq.size()>=param.min_read_size)&&(CountNs()<=param.max_ns)&&RunAlign(ref)) {
				StringAlign(ref, _str_align);
				n_aligned++;
				continue;
			}
			_pread->seq.erase(0,1);
			_pread->qual.erase(0,1);			
			while(1) {
				_pread->seq.erase(_pread->seq.size()-(param.trim_lowQ-50), param.trim_lowQ-50);
				_pread->qual.erase(_pread->qual.size()-(param.trim_lowQ-50), param.trim_lowQ-50);						
				if(_pread->seq.size()<param.min_read_size)
					break;
				if((CountNs()<=param.max_ns)&&RunAlign(ref)) {
					StringAlign(ref, _str_align);
					n_aligned++;
					break;
				}
			}
		}		
	}
}

void SingleAlign::DetectSnpSites(bool chain, Hit *hit, RefSeq &ref)
{
	int i,j;
	_sites.clear();
	if(!chain) {
		for(i=0; i<FIXELEMENT; i++)
			if(_a.a= (bseq[hit->z][i].a^ref.bfa[hit->chr].s[hit->loc/12+i].a)&reg[hit->z][i].a)
				for(j=0; j<12; j++)
					if((_a.a>>22-j*2)&0x3)
						_sites.push_back(i*12+j-hit->z);
	}
	else {
		for(i=0; i<FIXELEMENT; i++)
			if(_a.a= (cbseq[hit->z][i].a^ref.bfa[hit->chr].s[hit->loc/12+i].a)&creg[hit->z][i].a)
				for(j=0; j<12; j++)
					if((_a.a>>22-j*2)&0x3)
						_sites.push_back(i*12+j-hit->z);
	}
}

void SingleAlign::SetFlag(char c)
{
	_setname=c;
}

//output align hits
void SingleAlign::StringAlign(RefSeq &ref, string &os)
{
	Reverse_Seq();
	Reverse_Qual();
	int ii, jj, sum, j;
	//snp align:
	for(ii=0; ii<=param.max_snp_num; ii++) {
		if(0==_cur_n_hit[ii]+_cur_n_chit[ii]) //no hit
			continue;
		if(1==_cur_n_hit[ii]+_cur_n_chit[ii]) //unique hit
			if(ii==0) //exact align:
				if(_cur_n_hit[ii])
					s_OutHit(0, 1, ii, &hits[ii][0], 0, ref, os);
				else
					s_OutHit(1, 1, ii, &chits[ii][0], 0, ref, os);
			else
				if(_cur_n_hit[ii])
					s_OutHit(0, 1, ii, &hits[ii][0], 1, ref, os);
				else
					s_OutHit(1, 1, ii, &chits[ii][0], 1, ref, os);
		else if(1==param.report_repeat_hits) {  //randomly pick up one
			sum=(_cur_n_hit[ii]+_cur_n_chit[ii]<param.max_num_hits? _cur_n_hit[ii]+_cur_n_chit[ii] :param.max_num_hits);
			jj=rand()%sum;
			if(jj<_cur_n_hit[ii])
				s_OutHit(0, sum, ii, &hits[ii][jj], 1, ref, os);
			else
				s_OutHit(1, sum, ii, &chits[ii][jj-_cur_n_hit[ii]], 1, ref, os);
		}
		else if(2==param.report_repeat_hits) {   //output all repeat hits
			sum=(_cur_n_hit[ii]+_cur_n_chit[ii]<param.max_num_hits? _cur_n_hit[ii]+_cur_n_chit[ii] :param.max_num_hits);
			sort(hits[ii], hits[ii]+_cur_n_hit[ii], HitComp);	
			for(j=0; j<_cur_n_hit[ii]; j++) {
//				cout<<sum<<"  +  "<<(int)hits[ii][j].chr<<"  "<<hits[ii][j].loc<<"\n";
				s_OutHit(0, sum, ii, &hits[ii][j], 1, ref, os);
			}
			sort(chits[ii], chits[ii]+_cur_n_chit[ii], HitComp);
			for(j=0; j<_cur_n_chit[ii]; j++) {
//				cout<<sum<<"  -  "<<(int)chits[ii][j].chr<<"  "<<chits[ii][j].loc<<"\n";
				s_OutHit(1, sum, ii, &chits[ii][j], 1, ref, os);
			}
		}
		return;
	}
	//gap align:
	if(0==_cur_n_gaphit+_cur_n_cgaphit)
		return;
	if(1==_cur_n_gaphit+_cur_n_cgaphit) {
		if(_cur_n_gaphit)
			s_OutGapHit(0, 1, _gap_size, &gaphits[0], ref, os);
		else
			s_OutGapHit(1, 1, _gap_size, &cgaphits[0], ref, os);
	}
	else if(1==param.report_repeat_hits) {
		sum=(_cur_n_gaphit+_cur_n_cgaphit<param.max_num_hits? _cur_n_gaphit+_cur_n_cgaphit :param.max_num_hits);
		jj=rand()%sum;
		if(jj<_cur_n_gaphit)
			s_OutGapHit(0, sum, _gap_size, &gaphits[jj], ref, os);
		else
			s_OutGapHit(1, sum, _gap_size, &cgaphits[jj-_cur_n_gaphit], ref, os);
	}
	else if(2==param.report_repeat_hits) {
		sum=(_cur_n_gaphit+_cur_n_cgaphit<param.max_num_hits? _cur_n_gaphit+_cur_n_cgaphit :param.max_num_hits);
		sort(gaphits, gaphits+_cur_n_gaphit, HitComp);
		for(j=0; j<_cur_n_gaphit; j++)
			s_OutGapHit(0, sum, _gap_size, &gaphits[j], ref, os);
		sort(cgaphits, cgaphits+_cur_n_cgaphit, HitComp);
		for(j=0; j<_cur_n_cgaphit; j++)
			s_OutGapHit(1, sum, _gap_size, &cgaphits[j], ref, os);
	}
	return;
}

//write output according to types of hits
/* n: # of hits; chain: 0+/1-; flag: class of read; sig: 1, detect snp sites, 0, not */
void SingleAlign::s_OutHit(int chain, size_t n, bit8_t nsnps, Hit *hit, bool sig, RefSeq &ref, string &os)
{
	if(param.output_id)
		sprintf(_ch, "%s\t", _pread->name.c_str());
	else
		sprintf(_ch, "%u\t", _pread->index);
	os.append(_ch);
	if(!chain)
		sprintf(_ch, "%s\t%s\t", _pread->seq.c_str(), _pread->qual.c_str());
	else
		sprintf(_ch, "%s\t%s\t", _revseq.c_str(), _revqual.c_str());
	os.append(_ch);
	
	sprintf(_ch, "%d\t%c\t%d\t%c\t%s\t%u\t%d", n, _setname, _pread->seq.size(), chain_flag[chain], ref.title[hit->chr].name.c_str(), hit->loc+1, nsnps);
	os.append(_ch);
	if(sig && nsnps)
	{
		int i;
		DetectSnpSites(chain, hit, ref);
		if(!chain)
  		for(i=0; i!=_sites.size(); i++) {
  			sprintf(_ch, "\t%c->%d%c%d", nt_code[(ref.bfa[hit->chr].s[(hit->loc+_sites[i])/12].a>>((11-(hit->loc+_sites[i])%12)*2))&0x3], _sites[i], _pread->seq[_sites[i]], _pread->qual[_sites[i]]-param.zero_qual);
  			os.append(_ch);
  		}
  	else
  		for(i=0; i!=_sites.size(); i++) {
  			sprintf(_ch, "\t%c->%d%c%d", nt_code[(ref.bfa[hit->chr].s[(hit->loc+_sites[i])/12].a>>((11-(hit->loc+_sites[i])%12)*2))&0x3], _sites[i], _revseq[_sites[i]], _revqual[_sites[i]]-param.zero_qual);
  			os.append(_ch);
  		}  		
	}
	os.append("\n");
}

/* n: number of hits; g: gap size
*/
void SingleAlign::s_OutGapHit(int chain, size_t n, bit8_t g, Hit *hit, RefSeq &ref, string &os)
{
	if(param.output_id)
		sprintf(_ch, "%s\t", _pread->name.c_str());
	else
		sprintf(_ch, "%u\t", _pread->index);
	os.append(_ch);
	if(!chain)
		sprintf(_ch, "%s\t%s\t", _pread->seq.c_str(), _pread->qual.c_str());
	else
		sprintf(_ch, "%s\t%s\t", _revseq.c_str(), _revqual.c_str());
	os.append(_ch);
	
	//insertion on query
	if((hit->z <200)&&(hit->z>100))
		sprintf(_ch, "%d\t%c\t%d\t%c\t%s\t%u\t%d\t%d\n", n, _setname, _pread->seq.size(), chain_flag[chain], ref.title[hit->chr].name.c_str(), hit->loc+1, 100+g, hit->z-100);
	else if(hit->z>200)  //deletion on query
		sprintf(_ch, "%d\t%c\t%d\t%c\t%s\t%u\t%d\t%d\n", n, _setname, _pread->seq.size(), chain_flag[chain], ref.title[hit->chr].name.c_str(), hit->loc+1, 200+g, hit->z-200);
	else  //normal identical hit maybe
		sprintf(_ch, "%d\t%c\t%d\t%c\t%s\t%u\t%d\n", n, _setname, _pread->seq.size(), chain_flag[chain], ref.title[hit->chr].name.c_str(), hit->loc+1, 0);
	os.append(_ch);
}
