#include "param.h"

using namespace std;

Param::Param()
{
	num_procs=1;
	
	chains=0;

/*	
#ifdef DB_CHR  // seqs <256, length <4Gb
	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
#endif
#ifdef DB_CONTIG // seqs <65K, length <4Gb
	max_dbseq_size=0x100000; //1Mb
	append_dbseq_size=0x100000;  //1Mb
#endif
#ifdef DB_SHORT // seqs <4G, length <65K
	max_dbseq_size=0x10000; //65Kb
	append_dbseq_size=0x10000;  //65Kb
#endif
#ifdef DB_HUGE // seqs <4G, length <4G
	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
#endif
*/
	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
	
//	read_size=30;
	max_ns = 5;
	trim_lowQ=0;
	
	zero_qual= '@';
	qual_threshold= 20;
	default_qual=40;
	
	min_insert= 400;
	max_insert= 600;
	optimize_output_SV=1;
	
	seed_size= 10;
	half_seed_size= seed_size>>1;	
	half_seed_bits= (1<<(half_seed_size*2))-1;
	seed_bits=(1<<(seed_size*2))-1;
	min_read_size=half_seed_size*4+3;
	
	max_snp_num = 2;
	max_gap_size = 0;
	gap_edge = 5;
	max_num_hits = MAXHITS;
	
	//for mRNA tag alignment
	SetMrnaTag(-1);

	//for miRNA alignment
	adapter.clear();
	admis=0;
	mirna_min=17;
	mirna_max=26;
		
	report_repeat_hits = 1;
	output_id=1;

	useful_nt="ACGTacgt";
	nx_nt="NXnx";
	
	BuildMismatchTable();
};

void Param::SetSeedSize(int n)
{
	seed_size=(n/2)*2;
	half_seed_size= seed_size>>1;	
	half_seed_bits= (1<<(half_seed_size*2))-1;
	seed_bits=(1<<(seed_size*2))-1;
	min_read_size=half_seed_size*4+3;	
}

void Param::SetMrnaTag(int n)
{
	tag_type=n;
	if(-1==n)
		tag_seq.clear();
	else if(0==n) {
		tag_seq="GATC";
		tag_remain=16;
	}
	else if(1==n) {
		tag_seq="CATG";
		tag_remain=17;
	}
	else {
		cerr<<"fatal error: unrecognizable tag_type: "<<n<<endl;
		exit(1);
	}
}

void Param::BuildMismatchTable()
{
	int n;
	for(bit32_t i=0; i<0x1000000; i++) {
		n=0;
		for(int j=0; j<12; j++)
			if((i>>j*2)&0x3)
				n++;
		num_mismatch[i]=n;
	}
}
/*
Aa=0
Cc=1
Gg=2
Tt=3  */
const bit8_t nv=0;  //convert unknown char as 'A'
bit8_t alphabet[256] =
{
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'A' */ 
 0,nv, 1,nv,nv,nv, 2,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'P' */
nv,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'a' */ 
 0,nv, 1,nv,nv,nv, 2,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'p' */
nv,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv
};

bit8_t reg_alphabet[256] =
{
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'A' */ 
 3,nv, 3,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'P' */
nv,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'a' */ 
 3,nv, 3,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'p' */
nv,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv
};

const bit8_t nv2=3;  //'T'
bit8_t rev_alphabet[256] =
{
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2, /* next is 'A' */ 
3,nv2,2,nv2,nv2,nv2,1,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2, /* next is 'P' */
nv2,nv2,nv2,nv2,0,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2, /* next is 'a' */ 
3,nv2,2,nv2,nv2,nv2,1,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2, /* next is 'p' */
nv2,nv2,nv2,nv2,0,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,
nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2,nv2
};

const bit8_t nv3='N'; //set any unknown char as 'N'
char rev_char[256] =
{
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3, /* next is 'A' */ 
'T',nv3,'G',nv3,nv3,nv3,'C',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'P' */
nv3,nv3,nv3,nv3,'A',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3, /* next is 'a' */ 
't',nv3,'g',nv3,nv3,nv3,'c',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'p' */
nv3,nv3,nv3,nv3,'a',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3
};

const bit8_t nv4=4;
bit8_t n_char[256]=
{
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4, /* next is 'A' */ 
0,nv4,1,nv4,nv4,nv4,2,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4, /* next is 'P' */
nv4,nv4,nv4,nv4,3,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4, /* next is 'a' */ 
0,nv4,1,nv4,nv4,nv4,2,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4, /* next is 'p' */
nv4,nv4,nv4,nv4,3,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4
};
bit8_t nrev_char[256]=
{
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4, /* next is 'A' */ 
3,nv4,2,nv4,nv4,nv4,1,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4, /* next is 'P' */
nv4,nv4,nv4,nv4,0,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4, /* next is 'a' */ 
3,nv4,2,nv4,nv4,nv4,1,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4, /* next is 'p' */
nv4,nv4,nv4,nv4,0,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,
nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4,nv4
};
char chain_flag[2] =
{
'+', '-'
};

char nt_code[4] =
{
	'A', 'C', 'G', 'T'
};

char revnt_code[4] =
{
	'T', 'G', 'C', 'A'
};
