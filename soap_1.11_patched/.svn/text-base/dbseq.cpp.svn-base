#include "dbseq.h"

using namespace std;

extern Param param;
extern bit8_t alphabet[];

/************/
RefSeq::RefSeq()
{
	total_kmers=0;
}
ref_loc_t RefSeq::LoadNextSeq(ifstream &fin)
{
	char ch[1000];
	char c;
	string s;
	fin>>c;
	if(fin.eof())
		return 0;
	string::iterator z=_seq.begin();
	_length=0;
	//get name
	fin>>_name;
//	cout<<"name: "<<_name<<endl;
	fin.getline(ch, 1000);
	//get seq
	while(!fin.eof()) {
		fin>>c;
		if(fin.eof())
			break;
		fin.unget();
		if(c == '>')
			break;
		fin>>s;
		if(_length+s.size()>=param.max_dbseq_size) {
			param.max_dbseq_size+=param.append_dbseq_size;
			_seq.resize(param.max_dbseq_size);
			z=_seq.begin()+_length;
//			cout<<"_seq size: "<<param.max_dbseq_size<<endl;
		}
		copy(s.begin(), s.end(), z);
		z+=s.size();
		_length+=s.size();
	}
//	cout<<"got: "<<_length<<endl;
	return _length;
}

void RefSeq::BinSeq(OneBfa &a)
{
	a.n=(_length+11)/12+2;   //12bp, bit24 for each element. put 2 extra elements at the 3'end to invoid overflow
	int t=a.n*12-_length;
	if(t) {
		string ts(t, 'N');
		if(_seq.size()<_length+t)
			_seq.resize(_length+t);
		copy(ts.begin(), ts.end(), _seq.begin()+_length);
	}
	a.s = new bit24_t[a.n];
	bit32_t i=0;
	string::iterator p=_seq.begin();
	string::iterator tmp;
	for(; i<a.n; i++,p+=12) {
		a.s[i].a=0;
		for(bit32_t j=0; j<12; j++) {
			a.s[i].a<<=2;
			a.s[i].a|=alphabet[*(p+j)];
		}
	}
}

void RefSeq::UnmaskRegion()
{
	Block b;
	b.id=_count;
	b.begin=b.end=0;
//	bit32_t total_size=0;
	while(b.end<_length) {
		b.begin=_seq.find_first_of(param.useful_nt, b.end);
		if(b.begin > _length)
			break;
		b.end=_seq.find_first_of(param.nx_nt, b.begin);
		b.end = (b.end<=_length? b.end : _length);
		if(b.end-b.begin <30)
			continue;
		if((!_blocks.empty()) && (b.id==_blocks[_blocks.size()-1].id) 
			&& (b.begin - _blocks[_blocks.size()-1].end <5))
			_blocks[_blocks.size()-1].end=b.end;
		else {
			_blocks.push_back(b);
		}
	}
}

void RefSeq::Run_ConvertBinseq(ifstream &fin)
{
	_seq.resize(param.max_dbseq_size);
	RefTitle r;
	_count=0;
	total_num=sum_length=0;
	while(LoadNextSeq(fin)) {
		r.name=_name;
		r.size=_length;
		title.push_back(r);
		
		OneBfa a;
		BinSeq(a);
		bfa.push_back(a);
		UnmaskRegion();
		_count++;
		total_num++;
		sum_length+=_length;
//		cout<<r.size<<endl;
	}
//	cout<<"total seq length: "<<sum_length<<endl;
	_seq.clear(); //free ram
	//
/*	bit32_t a[25];
	bit32_t sum=0;
	for(int i=0; i<25; i++)
		a[i]=0;
	for(vector<Block>::iterator p=_blocks.begin(); p!=_blocks.end(); p++) {
		a[p->id]+=p->end-p->begin;
	}
	for(int i=0; i<25; i++) {
		cout<<i<<"  "<<a[i]<<endl;
		sum+=a[i];
	}
	cout<<"sum: "<<sum<<endl;*/
}

//continues seed
inline bit32_t RefSeq::s_MakeSeed_1(bit24_t *_m, int _a)
{
	bit32_t _seed;
	_seed= _a>=24?(_m->a>>(_a-24))&param.seed_bits:((_m->a<<(24-_a))|((_m+1)->a>>_a))&param.seed_bits;
	return _seed;
}
//spaced seed with one interval
inline bit32_t RefSeq::s_MakeSeed_2(bit24_t *_m, int _a, int _b, int _c, int _d)
{
	bit32_t _seed;
	if(_a>=24)
		if(_c>=24)
			_seed=(((_m->a>>(_a-24))&param.half_seed_bits)<<_d)|(((_m+_b)->a>>(_c-24))&param.half_seed_bits);
		else
			_seed=(((_m->a>>(_a-24))&param.half_seed_bits)<<_d)|(((_m+_b)->a<<(24-_c)|(_m+_b+1)->a>>_c)&param.half_seed_bits);
	else
		if(_c>=24)
			_seed=(((_m->a<<(24-_a)|(_m+1)->a>>_a)&param.half_seed_bits)<<_d)|(((_m+_b)->a>>(_c-24))&param.half_seed_bits);
		else
			_seed=(((_m->a<<(24-_a)|(_m+1)->a>>_a)&param.half_seed_bits)<<_d)|(((_m+_b)->a<<(24-_c)|(_m+_b+1)->a>>_c)&param.half_seed_bits);
	return _seed;
}
void RefSeq::InitialIndex()
{
	total_kmers=1<<param.seed_size*2;
	cout<<"total_kmers: "<<total_kmers<<endl;
	index= new KmerLoc[total_kmers];
	KmerLoc *p=index;
	bit32_t i=0;
	for(; i<total_kmers; p++,i++) {
		p->n1=p->n2=p->n3=0;
	}	
}

void RefSeq::t_CalKmerFreq_ab()
{
	bit24_t *_m;
	bit32_t i,e;
	int _a=48-param.seed_size*2;
	int _b=48-param.seed_size*2-8;
	int _c=48-param.seed_size*2-16;
	for(vector<Block>::iterator p=_blocks.begin(); p!=_blocks.end(); p++) {
		i=(p->begin+11)/12;  //
		e=(p->end-param.half_seed_size*2)/12;
		for(; i<=e; i++) {
			_m=bfa[p->id].s+i;
			index[s_MakeSeed_1(_m, _a)].n1++;
			index[s_MakeSeed_1(_m, _b)].n1++;
			index[s_MakeSeed_1(_m, _c)].n1++;
		}
	}
}
void RefSeq::t_CalKmerFreq_ac()
{
	bit24_t *_m;
	bit32_t i,e;
	int _a1=48-param.half_seed_size*2;
	int _b1=(param.half_seed_size*2)/12;
	int _c1=48-(param.half_seed_size+param.half_seed_size*2%12)*2;
	int _d1=param.half_seed_size*2;
	
	int _a2=48-param.half_seed_size*2-8;
	int _b2=(param.half_seed_size*2+4)/12;
	int _c2=48-(param.half_seed_size+(param.half_seed_size*2+4)%12)*2;
	int _d2=param.half_seed_size*2;
	
	int _a3=48-param.half_seed_size*2-16;
	int _b3=(param.half_seed_size*2+8)/12;
	int _c3=48-(param.half_seed_size+(param.half_seed_size*2+8)%12)*2;
	int _d3=param.half_seed_size*2;		
	for(vector<Block>::iterator p=_blocks.begin(); p!=_blocks.end(); p++) {
		i=(p->begin+11)/12;  //
		e=(p->end-param.half_seed_size*3)/12;
		for(; i<=e; i++) {
			_m=bfa[p->id].s+i;
			index[s_MakeSeed_2(_m, _a1, _b1, _c1, _d1)].n2++;
			index[s_MakeSeed_2(_m, _a2, _b2, _c2, _d2)].n2++;
			index[s_MakeSeed_2(_m, _a3, _b3, _c3, _d3)].n2++;						
		}
	}	
}
void RefSeq::t_CalKmerFreq_ad()
{
	bit24_t *_m;
	bit32_t i,e;
	int _a1=48-param.half_seed_size*2;
	int _b1=(param.half_seed_size*3)/12;
	int _c1=48-(param.half_seed_size+param.half_seed_size*3%12)*2;
	int _d1=param.half_seed_size*2;

	int _a2=48-param.half_seed_size*2-8;
	int _b2=(param.half_seed_size*3+4)/12;
	int _c2=48-(param.half_seed_size+(param.half_seed_size*3+4)%12)*2;
	int _d2=param.half_seed_size*2;
	
	int _a3=48-param.half_seed_size*2-16;
	int _b3=(param.half_seed_size*3+8)/12;
	int _c3=48-(param.half_seed_size+(param.half_seed_size*3+8)%12)*2;
	int _d3=param.half_seed_size*2;	
	for(vector<Block>::iterator p=_blocks.begin(); p!=_blocks.end(); p++) {
		i=(p->begin+11)/12;  //
		e=(p->end-param.half_seed_size*4)/12;
		for(; i<=e; i++) {
			_m=bfa[p->id].s+i;
			index[s_MakeSeed_2(_m, _a1, _b1, _c1, _d1)].n3++;
			index[s_MakeSeed_2(_m, _a2, _b2, _c2, _d2)].n3++;
			index[s_MakeSeed_2(_m, _a3, _b3, _c3, _d3)].n3++;						
		}
	}		
}
void RefSeq::AllocIndex()
{
	KmerLoc *v;
	bit32_t j;
	for(v=index, j=0; j<total_kmers; v++,j++)
		if(v->n1>0) {
			v->id1= new ref_id_t[v->n1];
			v->loc1= new ref_loc_t[v->n1];
			v->n1=0;
		}
	for(v=index, j=0; j<total_kmers; v++,j++)
		if(v->n2>0) {
			v->id2= new ref_id_t[v->n2];
			v->loc2= new ref_loc_t[v->n2];
			v->n2=0;
		}
	for(v=index, j=0; j<total_kmers; v++,j++)
		if(v->n3>0) {
			v->id3= new ref_id_t[v->n3];
			v->loc3= new ref_loc_t[v->n3];
			v->n3=0;
		}
}
void RefSeq::ReleaseIndex()
{
	for(bit32_t j=0; j<total_kmers; j++) {
		delete index[j].id1;
		delete index[j].loc1;
		delete index[j].id2;
		delete index[j].loc2;
		delete index[j].id3;
		delete index[j].loc3;
	}
	delete index;
}
void RefSeq::t_CreateIndex_ab()
{
	bit24_t *_m;
	KmerLoc *z;
	bit32_t i,e;
	int _a=48-param.seed_size*2;
	int _b=48-param.seed_size*2-8;
	int _c=48-param.seed_size*2-16;
	for(vector<Block>::iterator p=_blocks.begin(); p!=_blocks.end(); p++) {
		i=(p->begin+11)/12;  //
		e=(p->end-param.half_seed_size*2)/12;
		for(; i<=e; i++) {
			_m=bfa[p->id].s+i;
			z=index+s_MakeSeed_1(_m, _a);
			z->loc1[z->n1]=i*3;
			z->id1[z->n1]=p->id;
			z->n1++;
			z=index+s_MakeSeed_1(_m, _b);
			z->loc1[z->n1]=i*3+1;
			z->id1[z->n1]=p->id;
			z->n1++;
			z=index+s_MakeSeed_1(_m, _c);
			z->loc1[z->n1]=i*3+2;
			z->id1[z->n1]=p->id;
			z->n1++;						
		}
	}		
}
void RefSeq::t_CreateIndex_ac()
{
	bit24_t *_m;	
	KmerLoc *z;
	bit32_t i,e;
	int _a1=48-param.half_seed_size*2;
	int _b1=(param.half_seed_size*2)/12;
	int _c1=48-(param.half_seed_size+param.half_seed_size*2%12)*2;
	int _d1=param.half_seed_size*2;
	
	int _a2=48-param.half_seed_size*2-8;
	int _b2=(param.half_seed_size*2+4)/12;
	int _c2=48-(param.half_seed_size+(param.half_seed_size*2+4)%12)*2;
	int _d2=param.half_seed_size*2;
	
	int _a3=48-param.half_seed_size*2-16;
	int _b3=(param.half_seed_size*2+8)/12;
	int _c3=48-(param.half_seed_size+(param.half_seed_size*2+8)%12)*2;
	int _d3=param.half_seed_size*2;	
	for(vector<Block>::iterator p=_blocks.begin(); p!=_blocks.end(); p++) {
		i=(p->begin+11)/12;  //
		e=(p->end-param.half_seed_size*3)/12;
		for(; i<=e; i++) {
			_m=bfa[p->id].s+i;
			z=index+s_MakeSeed_2(_m, _a1, _b1, _c1, _d1);
			z->loc2[z->n2]=i*3;
			z->id2[z->n2]=p->id;
			z->n2++;
			z=index+s_MakeSeed_2(_m, _a2, _b2, _c2, _d2);
			z->loc2[z->n2]=i*3+1;
			z->id2[z->n2]=p->id;
			z->n2++;
			z=index+s_MakeSeed_2(_m, _a3, _b3, _c3, _d3);
			z->loc2[z->n2]=i*3+2;
			z->id2[z->n2]=p->id;
			z->n2++;						
		}
	}		
}
void RefSeq::t_CreateIndex_ad()
{
	bit24_t *_m;
	KmerLoc *z;
	bit32_t i,e;
	int _a1=48-param.half_seed_size*2;
	int _b1=(param.half_seed_size*3)/12;
	int _c1=48-(param.half_seed_size+param.half_seed_size*3%12)*2;
	int _d1=param.half_seed_size*2;

	int _a2=48-param.half_seed_size*2-8;
	int _b2=(param.half_seed_size*3+4)/12;
	int _c2=48-(param.half_seed_size+(param.half_seed_size*3+4)%12)*2;
	int _d2=param.half_seed_size*2;
	
	int _a3=48-param.half_seed_size*2-16;
	int _b3=(param.half_seed_size*3+8)/12;
	int _c3=48-(param.half_seed_size+(param.half_seed_size*3+8)%12)*2;
	int _d3=param.half_seed_size*2;	
	for(vector<Block>::iterator p=_blocks.begin(); p!=_blocks.end(); p++) {
		i=(p->begin+11)/12;  //
		e=(p->end-param.half_seed_size*4)/12;
		for(; i<=e; i++) {
			_m=bfa[p->id].s+i;
			z=index+s_MakeSeed_2(_m, _a1, _b1, _c1, _d1);
			z->loc3[z->n3]=i*3;
			z->id3[z->n3]=p->id;
			z->n3++;
			z=index+s_MakeSeed_2(_m, _a2, _b2, _c2, _d2);
			z->loc3[z->n3]=i*3+1;
			z->id3[z->n3]=p->id;
			z->n3++;
			z=index+s_MakeSeed_2(_m, _a3, _b3, _c3, _d3);
			z->loc3[z->n3]=i*3+2;
			z->id3[z->n3]=p->id;
			z->n3++;						
		}
	}
}
void RefSeq::CreateIndex()
{
	InitialIndex();
	t_CalKmerFreq_ab();
	t_CalKmerFreq_ac();
	t_CalKmerFreq_ad();	
	AllocIndex();
	t_CreateIndex_ab();	
	t_CreateIndex_ac();
	t_CreateIndex_ad();
}
