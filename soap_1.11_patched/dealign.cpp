#include "dealign.h"

using namespace std;

extern bit8_t reg_alphabet[];
extern bit8_t alphabet[];
extern bit8_t n_char[];
extern bit8_t nrev_char[];

void Dealign::IniCount(char *ref, vector<vector<bit16_t> > &c)
{
	ifstream fin(ref);
	if(!fin) {
		cerr<<"fatal error: failed to open ref "<<ref<<endl;
		exit(1);
	}
	_reftitle.clear();
	c.clear();
	vector<bit16_t> a;
	string s;
	char ch[1000];
	int i;
	while(!fin.eof()) {
		fin>>s;
		fin.getline(ch, 1000);
		if(fin.eof()) {
			if(!a.empty())
				c.push_back(a);
			break;
		}
		if(s[0] == '>') {
			s.erase(0,1);
			_reftitle.push_back(s);
			if(!a.empty())
				c.push_back(a);
			a.clear();
		}
		else {
			for(i=0; i!=s.size(); i++) {
				if(reg_alphabet[s[i]])
					a.push_back(1);
				else
					a.push_back(0);
			}
		}
	}
	fin.close();
	_titleorder.clear();
	vector<string>::iterator p=_reftitle.begin();
	for(i=0; p!=_reftitle.end(); p++,i++) {
		_titleorder[*p]=i;
	}
	//
	return;
}
void Dealign::OrderTitle()
{
	vector<string>::iterator p=_reftitle.begin();
	int i=0;
	for(; p!=_reftitle.end(); p++, i++) {
		_titleorder[*p]=i;
	}
	return;
}
vector<bit8_t> Dealign::GCcontent(string s, int win, int skip)
{
	vector<bit8_t> a;
	int gc, all;
	int i,j;
	for(i=0; i<s.size()-win; i+=skip) {
		all=gc=0;
		for(j=0; j!=win; j++) {
			if(reg_alphabet[s[i+j]])
				all++;
			if((alphabet[s[i+j]]==1) ||(alphabet[s[i+j]]==2))
				gc++;
		}
		if(all*2<win)   //too few meaningful sites
			a.push_back(0);
		else
			a.push_back((100*gc+all/2)/all);
	}
	return a;
}
void Dealign::CalGC(char *ref, int win, int slip)
{
	ifstream fin(ref);
	if(!fin) {
		cerr<<"fatal error: failed to open ref "<<ref<<endl;
		exit(1);
	}
	_gc.clear();
	string s;
	string ts;
	char ch[1000];
	while(!fin.eof()) {
		fin>>ts;
		fin.getline(ch, 1000);
		if(fin.eof()) {
			if(s.size())
				_gc.push_back(GCcontent(s,win,slip));
			break;
		}
		if(ts[0]=='>') {
			if(s.size())
				_gc.push_back(GCcontent(s,win,slip));
			s.clear();
		}
		else {
			s+=ts;
		}
	}
	fin.close();
}
/* flag==0, include all hits;
   flag==1, include only none repeat hits;
*/
void Dealign::CountHeadFreq(char *align_file, int flag)
{
	ifstream fin(align_file);
	if(!fin) {
		cerr<<"fatal error: failed to open align "<<align_file<<endl;
		exit(1);
	}	
	char ch[1000];
	int len, loc;
	char chain;
	AlignInfo a;
	string c_str;
	while(!fin.eof()) {
		fin>>a.id;
		if(fin.eof())
			break;
		fin>>a.seq>>a.qual>>a.nhits>>a.flag>>a.len>>a.chain>>c_str>>a.loc;
		fin.getline(a.note,NOTE_LEN);
		if(flag &&(a.nhits>1))
			continue;
		a.chr=_titleorder[c_str];
		if('+'==a.chain)
			_head[a.chr][a.loc]==max_depth? 0:_head[a.chr][a.loc]++;
		else
			_head[a.chr][a.loc+a.len-1]==max_depth? 0:_head[a.chr][a.loc+a.len-1]++;
	}
	fin.close();
}
/* flag==0, include all hits;
   flag==1, include only none repeat hits;
*/
void Dealign::CountDepth(char *align_file, int flag)
{
	ifstream fin(align_file);
	if(!fin) {
		cerr<<"fatal error: failed to open align "<<align_file<<endl;
		exit(1);
	}	
	char ch[1000];
	int len, loc;
	char chain;
	AlignInfo a;
	string c_str;
	int i;
	while(!fin.eof()) {
		fin>>a.id;
		if(fin.eof())
			break;
		fin>>a.seq>>a.qual>>a.nhits>>a.flag>>a.len>>a.chain>>c_str>>a.loc;
		fin.getline(a.note,NOTE_LEN);
		if(flag &&(a.nhits>1))
			continue;
		a.chr=_titleorder[c_str];
		for(i=0; i!=a.len; i++)
				_depth[a.chr][a.loc+i]==max_depth? 0:_depth[a.chr][a.loc+i]++;
	}
	fin.close();
}
void Dealign::OutHist(char *out_file, vector<vector<bit16_t> > &c, int win, int slip, float zoom)
{
	ofstream fout(out_file);
	if(!fout) {
		cerr<<"fatal error: failed to open output file "<<out_file<<endl;
		exit(1);
	}
	int i,j;
	vector<vector<bit16_t> >::iterator p;
	bit64_t freq[max_depth+1];
	for(i=0; i<max_depth; i++)
		freq[i]=0;
	bit64_t sum, num, zsum;
	if(win >1) {
		for(p=c.begin(); p!=c.end(); p++) {
			for(i=0; i<p->size()-win; i+=slip) {
				sum=num=0;
				for(j=0; j<win; j++) {
					if((*p)[i+j]!=0) {   //exclude 'N's on reference
						num++;
						sum+=(*p)[i+j];
					}
				}
				if(num>=win*0.5) {
					zsum=int(zoom*(sum-num)/num+0.5);
					freq[zsum>max_depth? max_depth: zsum]++; //now is real depth, by zoomed
				}
			}
		}
	}
	else {
		for(p=c.begin(); p!=c.end(); p++) {
			for(i=0; i!=p->size(); i+=slip) {
				if((*p)[i]!=0)
					freq[(*p)[i]-1]++; 
			}
		}
	}
	//
	sum=0;
	for(i=0; i!=max_depth; i++)
		sum+=freq[i];
	int max_d;
	for(i=max_depth; i>=0; i--)
		if(freq[i]) {
			max_d=i;
			break;
		}
	for(i=0; i<=max_d; i++)
		fout<<i<<"\t"<<freq[i]<<"\t"<<freq[i]/double(sum/100)<<"%\n";
	fout.close();
}
void Dealign::OutDistri(char *out_file, vector<vector<bit16_t> > &c, int win, int slip, float zoom)
{
	ofstream fout(out_file);
	if(!fout) {
		cerr<<"fatal error: failed to open output file "<<out_file<<endl;
		exit(1);
	}
	bit32_t i,j,k;
	vector<vector<bit16_t> >::iterator p;
	if(win >1) {
		int sum, num;
		for(p=c.begin(), k=0; p!=c.end(); p++,k++) {
			fout<<">"<<_reftitle[k]<<"  "<<(p->size()-win)/slip+1<<endl;
			for(i=0; i<p->size()-win; i+=slip) {
				sum=num=0;
				for(j=0; j<win; j++) {
					if((*p)[i+j]!=0) {   //exclude 'N's on reference
						num++;
						sum+=(*p)[i+j];
					}
				}
				if(num>=win*0.5)
					fout<<int(zoom*(sum-num)/num+0.5)<<" ";  //real depth*zoom
				else
					fout<<-1<<" ";
			}
			fout<<"\n";
		}
	}
	else {
		for(p=c.begin(), k=0; p!=c.end(); p++,k++) {
			fout<<">"<<_reftitle[k]<<"  "<<p->size()<<endl;
			for(i=0; i!=p->size(); i++) {
				if((*p)[i]>0)
					fout<<(*p)[i]-1<<" ";
				else
					fout<<-1<<" ";
			}
			fout<<endl;
		}
	}
	fout.close();
}
void Dealign::OutGCdot(char *out_file, vector<vector<bit16_t> > &c, int win, int slip, float zoom)
{
	ofstream fout(out_file);
	if(!fout) {
		cerr<<"fatal error: failed to open output file "<<out_file<<endl;
		exit(1);
	}
	bit32_t i,j,k,m;
	vector<vector<bit16_t> >::iterator p;
	if(win <=1) {
		cerr<<"fatal error: GC vs depth dot plot can not use 1-bp win\n";
		exit(1);
	}
	int sum, num;
	for(p=c.begin(), k=0; p!=c.end(); p++,k++) {
		for(i=0,m=0; i<p->size()-win; i+=slip,m++) {
			sum=num=0;
			for(j=0; j<win; j++) {
				if((*p)[i+j]!=0) {   //exclude 'N's on reference
					num++;
					sum+=(*p)[i+j];
				}
			}
			if(num>=win*0.5)
				fout<<i<<"  "<<(int)_gc[k][m]<<"  "<<int(zoom*(sum-num)/num+0.5)<<"\n";
			else
				fout<<i<<"  "<<(int)_gc[k][m]<<"  "<<-1<<"\n";
		}
	}
	fout.close();	
}
void Dealign::LoadChrOrder(char *chrorder_file)
{
	ifstream fin(chrorder_file);
	if(!fin) {
		cerr<<"fatal error: failed to open chr order file "<<chrorder_file<<endl;
		exit(1);
	}
	_reftitle.clear();
	_titleorder.clear();
	string c_str;
	char ch[1000];
	int i=0;
	while(!fin.eof()) {
		fin>>c_str;
		fin.getline(ch, 1000);
		if(fin.eof())
			break;
		_reftitle.push_back(c_str);
		_titleorder[c_str]=i++;
	}
	return;
}
void Dealign::LoadAlign(char *align_file)
{
	ifstream fin(align_file);
	if(!fin) {
		cerr<<"fatal error: failed to open align "<<align_file<<endl;
		exit(1);
	}		
	_align.clear();
	AlignInfo a;
	string c_str;
	while(!fin.eof()) {
		fin>>a.id;
		if(fin.eof())
			break;
		fin>>a.seq>>a.qual>>a.nhits>>a.flag>>a.len>>a.chain>>c_str>>a.loc;
		a.chr=_titleorder[c_str];
		fin.getline(a.note,NOTE_LEN);
		_align.push_back(a);
	}
}

bool less_chr_loc(AlignInfo a, AlignInfo b)
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
}

bool bigger_chr_loc(AlignInfo a, AlignInfo b)
{
	if(a.chr>b.chr)
		return 1;
	else if(a.chr==b.chr) {
		if(a.loc>b.loc)
			return 1;
		else
			return 0;
	}
	else
		return 0;
}
void Dealign::SortAlign(int flag)
{
	if(0==flag) //0->1
		sort(_align.begin(), _align.end(), less_chr_loc);
	else
		sort(_align.begin(), _align.end(), bigger_chr_loc);
	return;
}
void Dealign::ExportAlign(char *align_file)
{
	ofstream fout(align_file);
	if(!fout) {
		cerr<<"fatal error: failed to create file "<<align_file<<endl;
		exit(1);
	}
	for(vector<AlignInfo>::iterator p=_align.begin(); p!=_align.end(); p++) {
		fout<<p->id<<"\t"<<p->seq<<"\t"<<p->qual<<"\t"<<p->nhits<<"\t"<<p->flag<<"\t"<<p->len<<"\t"
		<<p->chain<<"\t"<<_reftitle[p->chr]<<"\t"<<p->loc<<p->note<<"\n";
	}
	fout.close();
	return;
}
void Dealign::MergeSorted(char *infiles, char *outfile)
{
}
void Dealign::IniQC()
{
	int i,j;
	for(i=0; i!=2; i++)
		for(j=0; j!=10; j++)
			_ur[i][j]=0;
	for(i=0; i!=SEQ_LEN; i++) {
		for(j=0; j!=5; j++)
			_ntfreq[i][j]=0;
		for(j=0; j!=6; j++)
			_mut[i][j]=0;
		for(j=0; j!=20; j++)
			_match_qual[i][j]=_mis_qual[i][j]=0;
	}
}

/*read_len: full length of original raw read,
  flag: 0, look at only full-length reads,
        1, look at all reads, didn't trim first bp
        2, look at all reads, trimmed first bp,
  zero_qual: zero quality, default='@' in solexa
*/
void MutDe(string mismatch, int &offset, char &allele, int &qvalue){
	unsigned int i = 3;
	if (isdigit(mismatch[i+1])){
		offset = (mismatch[i]-48) *10 + (mismatch[i+1]-48);
		i += 2;
		allele=mismatch[i++];
	}
	else if (!isdigit(mismatch[i+1])){
		offset = (mismatch[i]-48);
		i +=1;
		allele=mismatch[i++];
	}			
	if (isdigit(mismatch[i])&&isdigit(mismatch[i+1]))
		qvalue = (mismatch[i]-48)*10 + (mismatch[i+1] -48);
	else
		qvalue = (mismatch[i]-48);
}

void Dealign::QC(char *align_file, int read_len, int flag, char zero_qual)
{
	ifstream fin(align_file);
	if(!fin) {
		cerr<<"fatal error: failed to open align "<<align_file<<endl;
		exit(1);
	}
	string id,seq,qual,mismatch;
	int nhits,len,nmut,offset,qvalue;
	char chain,allele;
	char ch[1000];
	int i,j,k,t;
	
	if(0==flag) {
  	while(!fin.eof()) {
  		fin>>id;
  		if(fin.eof())
  			break;
  		fin>>seq>>qual>>nhits>>ch>>len>>chain>>ch>>ch>>nmut;
  		if(nmut>100) {   //gapped alignment, will skip in stastics for this momment
  			fin.getline(ch, 1000);
  			continue;
  		}
  		nhits==1?_ur[0][nmut]++ : _ur[1][nmut]++;
  		if(nhits>1) {     //repeat hits, ignore them for this momment too
  			fin.getline(ch, 1000);
  			continue;
  		}			
			if(len!=read_len)
				continue;
			if('+'==chain) {
				for(i=0; i!=seq.size(); i++) {
					_ntfreq[i][n_char[seq[i]]]++;
					qvalue=(qual[i]>=zero_qual? (qual[i]-zero_qual)/10:0);
					_match_qual[i][qvalue]++;
				}			
				if(nmut) {
					for(k=0; k!=nmut; k++) {
						fin>>mismatch;
						MutDe(mismatch, offset, allele, qvalue);
						_mut[offset][n_char[allele]]++;
						_mis_qual[offset][qvalue>=0? qvalue/10:0]++;
					}
				}
			}
			else {
				for(i=0, j=seq.size()-1; j>=0; i++,j--) {
					_ntfreq[i][nrev_char[seq[j]]]++;
					qvalue=(qual[j]>=zero_qual? (qual[j]-zero_qual)/10:0);
					_match_qual[i][qvalue]++;
				}
				if(nmut) {
					for(k=0; k!=nmut; k++) {
						fin>>mismatch;
						MutDe(mismatch, offset, allele, qvalue);
						_mut[seq.size()-1-offset][nrev_char[allele]]++;
						_mis_qual[seq.size()-1-offset][qvalue>=0? qvalue/10:0]++;
					}
				}
			}
			_mut[seq.size()-1][4]++;
		}
	}	
	else if(1==flag) {
  	while(!fin.eof()) {
  		fin>>id;
  		if(fin.eof())
  			break;
  		fin>>seq>>qual>>nhits>>ch>>len>>chain>>ch>>ch>>nmut;
  		if(nmut>100) {   //gapped alignment, will skip in stastics for this momment
  			fin.getline(ch, 1000);
  			continue;
  		}
  		nhits==1?_ur[0][nmut]++ : _ur[1][nmut]++;
  		if(nhits>1) {     //repeat hits, ignore them for this momment too
  			fin.getline(ch, 1000);
  			continue;
  		}				
			if('+'==chain) {
				for(i=0; i!=seq.size(); i++) {
					_ntfreq[i][n_char[seq[i]]]++;
					qvalue=(qual[i]>=zero_qual? (qual[i]-zero_qual)/10:0);
					_match_qual[i][qvalue]++;
				}
				if(nmut) {
					for(k=0; k!=nmut; k++) {
						fin>>mismatch;
						MutDe(mismatch, offset, allele, qvalue);
						_mut[offset][n_char[allele]]++;
						_mis_qual[offset][qvalue>=0? qvalue/10:0]++;
					}
				}
			}
			else {
				for(i=0, j=seq.size()-1; j>=0; i++,j--) {
					_ntfreq[i][nrev_char[seq[j]]]++;
					qvalue=(qual[j]>=zero_qual? (qual[j]-zero_qual)/10:0);
					_match_qual[i][qvalue]++;
				}
				if(nmut) {
					for(k=0; k!=nmut; k++) {
						fin>>mismatch;
						MutDe(mismatch, offset, allele, qvalue);
						_mut[seq.size()-1-offset][nrev_char[allele]]++;
						_mis_qual[seq.size()-1-offset][qvalue>=0? qvalue/10:0]++;
					}
				}
			}
			_mut[seq.size()-1][4]++;			
		}
	}
		
	else if(2==flag) {
  	while(!fin.eof()) {
  		fin>>id;
  		if(fin.eof())
  			break;
  		fin>>seq>>qual>>nhits>>ch>>len>>chain>>ch>>ch>>nmut;
  		if(nmut>100) {   //gapped alignment, will skip in stastics for this momment
  			fin.getline(ch, 1000);
  			continue;
  		}
  		nhits==1?_ur[0][nmut]++ : _ur[1][nmut]++;
  		if(nhits>1) {     //repeat hits, ignore them for this momment too
  			fin.getline(ch, 1000);
  			continue;
  		}			
			t=(seq.size()==read_len? 0:1);
			if('+'==chain) {
				for(i=0; i!=seq.size(); i++) {
					_ntfreq[i+t][n_char[seq[i]]]++;
					qvalue=(qual[i]>=zero_qual? (qual[i]-zero_qual)/10:0);
					_match_qual[i+t][qvalue]++;
				}
				if(nmut) {
					for(k=0; k!=nmut; k++) {
						fin>>mismatch;
						MutDe(mismatch, offset, allele, qvalue);
						_mut[offset+t][n_char[allele]]++;
						_mis_qual[offset+t][qvalue>=0? qvalue/10:0]++;
					}
				}
			}
			else {
				for(i=0, j=seq.size()-1; j>=0; i++,j--) {
					_ntfreq[i+t][nrev_char[seq[j]]]++;
					qvalue=(qual[j]>=zero_qual? (qual[j]-zero_qual)/10:0);
					_match_qual[i+t][qvalue]++;
				}
				if(nmut) {
					for(k=0; k!=nmut; k++) {
						fin>>mismatch;
						MutDe(mismatch, offset, allele, qvalue);
						_mut[seq.size()-1-offset+t][nrev_char[allele]]++;
						_mis_qual[seq.size()-1-offset+t][qvalue>=0? qvalue/10:0]++;
					}
				}
			}
			_mut[seq.size()-1+t][4]++;			
		}		
		fin.getline(ch,1000);
	}
	fin.close();
}
void Dealign::OutQC(char *out_file)
{
	FILE * fout=fopen(out_file, "w");
	if(fout==NULL) {		
		cerr<<"fatal error: failed to open "<<out_file<<endl;
		exit(1);
	}

	int i,j,k,n;
	bit32_t total_mapped=0;
	for(i=0; i!=2; i++)
		for(j=0; j!=3; j++)
			total_mapped+=_ur[i][j];
	fprintf(fout, "Total aligned reads: %10d\n",total_mapped);
	fprintf(fout, "Unique, 0 mutation:  %10d, %5.1f%\n",_ur[0][0],100*double(_ur[0][0])/total_mapped);
	fprintf(fout, "Unique, 1 mutation:  %10d, %5.1f%\n",_ur[0][1],100*double(_ur[0][1])/total_mapped);
	fprintf(fout, "Unique, 2 mutation:  %10d, %5.1f%\n",_ur[0][2],100*double(_ur[0][2])/total_mapped);
	fprintf(fout, "Repeat, 0 mutation:  %10d, %5.1f%\n",_ur[1][0],100*double(_ur[1][0])/total_mapped);
	fprintf(fout, "Repeat, 1 mutation:  %10d, %5.1f%\n",_ur[1][1],100*double(_ur[1][1])/total_mapped);
	fprintf(fout, "Repeat, 2 mutation:  %10d, %5.1f%\n",_ur[1][2],100*double(_ur[1][2])/total_mapped);
	
	fprintf(fout, "\n");
	fprintf(fout, "Note: following only looked at uniquely aligned hits\n\n");
	fprintf(fout, "                     GC content                          Mutation (*1000)             Quality values          Quality of mutation sites \n");
	fprintf(fout, "        _______________________________             ________________________   ___________________________   ___________________________\n");
	fprintf(fout, "Base    A    C    G    T    N   GC   AT    #reads   rate    A    C    G    T   0~9 10~19 20~29 30~39 40~49   0~9 10~19 20~29 30~39 40~49\n");
	for(n=SEQ_LEN-1; n>=0; n--) {
		if(_mut[n][4])
			break;
	}
	double nreads, total_mut;
	for(i=0; i<=n; i++) {
		nreads=(_ntfreq[i][0]+_ntfreq[i][1]+_ntfreq[i][2]+_ntfreq[i][3]+_ntfreq[i][4])/100;
		//GC
		fprintf(fout, "%4d%5.1f%5.1f%5.1f%5.1f%5.1f%5.1f%5.1f",
		 i+1,_ntfreq[i][0]/nreads,_ntfreq[i][1]/nreads,_ntfreq[i][2]/nreads,_ntfreq[i][3]/nreads,
		 _ntfreq[i][4]/nreads,(_ntfreq[i][1]+_ntfreq[i][2])/nreads,(_ntfreq[i][0]+_ntfreq[i][3])/nreads);
		//mutation bias
		nreads/=10;
		total_mut=_mut[i][0]+_mut[i][1]+_mut[i][2]+_mut[i][3];
		fprintf(fout, "%10d  %5.1f%5.1f%5.1f%5.1f%5.1f",_mut[i][4],total_mut/nreads,
		 _mut[i][0]/nreads,_mut[i][1]/nreads,_mut[i][2]/nreads,_mut[i][3]/nreads);	 
		//base pair quality		
		fprintf(fout, "%6d%6d%6d%6d%6d",int(_match_qual[i][0]/nreads),int(_match_qual[i][1]/nreads),
			int(_match_qual[i][2]/nreads),int(_match_qual[i][3]/nreads));
		//mismatch quality
		total_mut/=1000;
		fprintf(fout, "%6d%6d%6d%6d%6d",int(_mis_qual[i][0]/total_mut),int(_mis_qual[i][1]/total_mut),
			int(_mis_qual[i][2]/total_mut),int(_mis_qual[i][3]/total_mut));
		fprintf(fout, "\n");
	}
	//
	fclose(fout);
}
