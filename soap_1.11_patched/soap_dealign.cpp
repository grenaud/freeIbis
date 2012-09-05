#include<unistd.h>
#include<cstdio>
#include<iostream>
#include<ostream>
#include<fstream>
#include<string>
#include<vector>
#include "dealign.h"

using namespace std;

void usage(void)
{
	cout<<"Usage: soap_dealign	[options]\n"
		<<"      sortalign  <in.chrorder> <in.align> <out.align>\n\n"
		<<"      headhist   <win> <slip> <zoom> <0:all,1:rm repeat hits> <ref.fa> <outfile> <in.align 1> <in.align 2> ...\n"
		<<"                 Output format: depth, frequency, % of all chr locations\n\n"
		<<"      depthhist  <win> <slip> <zoom> <0:all,1:rm repeat hits> <ref.fa> <outfile> <in.align 1> <in.align 2> ...\n"
		<<"                 Output format: depth, frequency, % of all chr locations\n\n"
		<<"      headdis    <win> <slip> <zoom> <0:all,1:rm repeat hits> <ref.fa> <outfile> <in.align 1> <in.align 2> ...\n\n"
		<<"      depthdis   <win> <slip> <zoom> <0:all,1:rm repeat hits> <ref.fa> <outfile> <in.align 1> <in.align 2> ...\n\n"
		<<"      gc2depth   <win> <slip> <zoom> <0:all,1:rm repeat hits> <ref.fa> <outfile> <in.align 1> <in.align 2> ...\n"
		<<"                 Output format: location, GC, mean depth\n\n"
		<<"      QC         <read_len> <flag> <char_zero_quality> <outfile> <in.align 1> <in.align 2> ...\n"
		<<"                 Note: flag=0, look only full-length reads; flag=1, look all reads; flag=2, look all reads(first bp trimmed align)\n"
		<<"                       gapped hits and repeat hits will be excluded from the calculation for this momment\n";
		exit(1);
};

int main(int argc, char *argv[])
{
	//print usage
	if ((argc==1) ||(argc==2))
	{
		usage();
	}
	Dealign a;
	if(!strcmp("sortalign", argv[1])) {  //sort alignment according to location on reference
		a.LoadChrOrder(argv[2]);
		a.LoadAlign(argv[3]);
		a.SortAlign(0);
		a.ExportAlign(argv[4]);
	}
	else if(!strcmp("headhist", argv[1])) {  //histogram of read head distribution
		a.IniCount(argv[6], a._head);
		for(int i=8; i!=argc; i++)
			a.CountHeadFreq(argv[i], atoi(argv[5]));
		a.OutHist(argv[7], a._head, atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
	}
	else if(!strcmp("depthhist", argv[1])) {  //histogram of read depth distribution
		a.IniCount(argv[6], a._depth);
		for(int i=8; i!=argc; i++)
			a.CountDepth(argv[i], atoi(argv[5]));
		a.OutHist(argv[7], a._depth, atoi(argv[2]), atoi(argv[3]), atof(argv[4]));		
	}
	else if(!strcmp("headdis", argv[1])) { //distribution of read head along reference
		a.IniCount(argv[6], a._head);
		for(int i=8; i!=argc; i++)
			a.CountHeadFreq(argv[i], atoi(argv[5]));
		a.OutDistri(argv[7], a._head, atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
	}
	else if(!strcmp("depthdis", argv[1])) {  //distribution of read depth along reference
		a.IniCount(argv[6], a._depth);
		for(int i=8; i!=argc; i++)
			a.CountDepth(argv[i], atoi(argv[5]));
		a.OutDistri(argv[7], a._depth, atoi(argv[2]), atoi(argv[3]), atof(argv[4]));		
	}
	else if(!strcmp("gc2depth", argv[1])) {  //gc vs depth in a window
		a.IniCount(argv[6], a._depth);
		a.CalGC(argv[6], atoi(argv[2]), atoi(argv[3]));
		for(int i=8; i!=argc; i++)
			a.CountDepth(argv[i], atoi(argv[5]));
		a.OutGCdot(argv[7], a._depth, atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
	}
	else if(!strcmp("QC", argv[1])) {  //summary and QC check of alignment
		a.IniQC();
		for(int i=6; i!=argc; i++)
			a.QC(argv[i], atoi(argv[2]), atoi(argv[3]), argv[4][0]);
		a.OutQC(argv[5]);
	}
	return 0;
}
