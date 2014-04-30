/*
 * combinePhixFiles
 * Date: Apr-22-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <sstream>      

#include "func.h"

using namespace std;

typedef struct{
    unsigned int cov;
    unsigned int base[4];//for mismatches only
    // unsigned int mmC;
    // unsigned int mmG;
    // unsigned int mmT;
} coverInfo;

inline int base2int(const char c){
    if(c ==    'A')
	return 0;
    if(c ==    'C')
	return 1;
    if(c ==    'G')
	return 2;
    if(c ==    'T')
	return 3;

    cerr<<"combinePhixFiles base2int() Invalid base "<<c<<endl;
    exit(1);
}

int main (int argc, char *argv[]) {

    double  fractionMask=1.0;
    bool nomask=false;
    const string usage=string(argv[0])+"\t [reference for phix] [bcl2phix output 1] [bcl2phix output 2] ...\n"+
	"Options:\n"+
	"\tMandatory\n"+
	"\t\t-f [fraction coverage phix for mask (Default: no masking)]\n"+
	"\t\t-m [mask output]\n"+
	"\t\t-c [coverage output]\n"+
	"\t\t-o [phix training output]\n"+
	"\tOptional\n"+
	"\t\t-n Do not mask positions on the PhiX (Default: will mask)\n"+
	"\n";
    //CGATTCG

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<usage<<endl;
	
	return 1;
    }

    ofstream maskOut;
    ofstream covOut;
    ofstream seqOut;
    string maskOutS;
    string covOutS;
    string seqOutS;
    int argCount=1;

    for(argCount=1;argCount<(argc);argCount++){

	if(string(argv[argCount])[0] != '-')
	    break;

	if(string(argv[argCount]) == "-f" ){
	    fractionMask=atof(argv[argCount+1]);
	    argCount++;
	    continue;
	}

	if(string(argv[argCount]) == "-m" ){
	    maskOutS=string(argv[argCount+1]);
	    argCount++;
	    continue;
	}

	if(string(argv[argCount]) == "-n" ){
	    nomask=true;
	    continue;
	}

	if(string(argv[argCount]) == "-c" ){
	    covOutS=string(argv[argCount+1]);
	    argCount++;
	    continue;
	}

	if(string(argv[argCount]) == "-o" ){
	    seqOutS=string(argv[argCount+1]);
	    argCount++;
	    continue;
	}


    }

    string phixFileName      = string(argv[argCount++]);
    if(!nomask)
	maskOut.open(maskOutS.c_str(), ofstream::out);
    covOut.open(covOutS.c_str(),   ofstream::out);
    seqOut.open(seqOutS.c_str(),   ofstream::out);

    if(!nomask)
	if( !maskOut.good() ){
	    cerr<<"Cannot write to file "<<maskOutS<<endl;
	}
    
    if( !covOut.good() ){
	cerr<<"Cannot write to file "<<covOutS<<endl;
    }

    if( !seqOut.good() ){
	cerr<<"Cannot write to file "<<seqOutS<<endl;
    }




    ifstream phixFile;


    

    string phixgenome   = "";
    string line;
    string namePhixRef;

    phixFile.open(phixFileName.c_str(), ios::in);

    if (phixFile.is_open()){
	getline (phixFile,line);

	namePhixRef=line.substr(1);

	if(namePhixRef.find(" ") != string::npos){
	    namePhixRef=namePhixRef.substr(0,line.find(" "));
	}

	
	while ( getline (phixFile,line)){
	    transform(line.begin(), line.end(), line.begin(), ::toupper);
	    phixgenome+=line;
	}
	phixFile.close();
    }else{
	cerr << "Unable to open file "<<phixFileName<<endl;
	return 1;
    }

    string phixgenomeRC=reverseComplement(phixgenome);

    vector<coverInfo> vecCov;
    for(unsigned int i=0;i<phixgenome.size();i++){
	coverInfo toadd;
	toadd.cov=0;
	for(unsigned int j=0;j<4;j++)
	    toadd.base[j]=0;

	vecCov.push_back(toadd);
    }


    for(int i=argCount;i<argc;i++){
	ifstream myFile;
	string filename = string(argv[i]);

	myFile.open(filename.c_str(), ios::in);

	if (myFile.is_open()){
	    while ( getline (myFile,line)){
		//cout<<line<<endl;
		stringstream   ss (line);
		string temp;
		// string temp1;

		string seq;
		string poss;
		int pos;

		for(int k=0;k<6;k++){
		    getline(ss, temp, '\t');
		    // temp1+=temp+"\t";
		}

		getline(ss, seq, '\t');
		getline(ss, poss, '\t');
		pos = atoi(poss.c_str());
		// seqOut<<temp1;
		
		if(pos < 0 ){ //- strand
		    // seqOut<<phixgenomeRC.substr(-1*pos,seq.size())<<endl;
		     // cout<<"2"<<phixgenome.substr(phixgenome.size()+pos-int(seq.size()),seq.size())<<endl;
		    //cout<<"seq "<<seq<<endl;
		    unsigned int posSeq=seq.size()-1;
		    pos++;//bring back to zero based
		    for(unsigned int posPhiX=(phixgenome.size()+pos-int(seq.size()));posPhiX<((phixgenome.size()+pos-int(seq.size()))+seq.size());posPhiX++,posSeq--){

			char b=complement(seq[posSeq]);
			//cout<<phixgenome[posPhiX]<<"\t"<<b<<endl;
			vecCov[posPhiX].cov++;
			if(phixgenome[posPhiX] != b &&
			   b                   != 'N'){
			    vecCov[posPhiX].base[ base2int( b) ]++;
			}
		    }

		}else{ //+strand
		    // seqOut<<phixgenome.substr(pos,seq.size())<<endl;
		    unsigned int posSeq=0;
		    pos--;//bring back to zero based
		    for(unsigned int posPhiX=pos;posPhiX<(pos+seq.size());posPhiX++,posSeq++){
			vecCov[posPhiX].cov++;
			//cout<<phixgenome[posPhiX]<<"\t"<< seq[posSeq]<<endl;
			if(phixgenome[posPhiX] != seq[posSeq] &&
			   seq[posSeq]         != 'N'){
			    vecCov[posPhiX].base[ base2int(seq[posSeq]) ]++;
			}
		    }

		}
		// }
	    }
	    myFile.close();
	}else{
	    cerr << "Unable to open file "<<filename<<endl;
	    return 1;
	}
    }

    vector<unsigned int> vecMaskPos;
    for(unsigned int i=0;i<phixgenome.size();i++){
	covOut<<namePhixRef<<"\t"<<(i)<<"\t"<<vecCov[i].cov<<"\t"<<vecCov[i].base[0]<<"\t"<<vecCov[i].base[1]<<"\t"<<vecCov[i].base[2]<<"\t"<<vecCov[i].base[3]<<endl;	
	if(vecCov[i].cov == 0)
	    continue;
	//cout<<i<<"\t"<<double(vecCov[i].base[0])/double(vecCov[i].cov)<<endl;
	if( (double(vecCov[i].base[0])/double(vecCov[i].cov)  > fractionMask) ||
	    (double(vecCov[i].base[1])/double(vecCov[i].cov)  > fractionMask) ||
	    (double(vecCov[i].base[2])/double(vecCov[i].cov)  > fractionMask) ||
	    (double(vecCov[i].base[3])/double(vecCov[i].cov)  > fractionMask) ){
	    
	    vecMaskPos.push_back(i);
	}
	    
    }


    //BEGIN MASKING
    if(!nomask){
	for(unsigned int i=0;i<vecMaskPos.size();i++){
	    maskOut<<namePhixRef<<"\t"<<vecMaskPos[i]<<endl;
	    phixgenome[vecMaskPos[i]] = 'N';
	    phixgenomeRC[phixgenome.size()-1-vecMaskPos[i]] = 'N';
	    
	}
	
	maskOut.close();
    }

    covOut.close();


    for(int i=argCount;i<argc;i++){
	ifstream myFile;
	string filename = string(argv[i]);

	myFile.open(filename.c_str(), ios::in);

	if (myFile.is_open()){
	    while ( getline (myFile,line)){
		// cout<<line<<endl;
		stringstream   ss (line);
		string temp;
		string temp1;

		string seq;
		string poss;
		int pos;

		for(int k=0;k<6;k++){
		    getline(ss, temp, '\t');
		    temp1+=temp+"\t";
		}

		getline(ss, seq, '\t');
		getline(ss, poss, '\t');
		pos = atoi(poss.c_str());
		seqOut<<temp1;
		
		if(pos < 0 ){
		    
		     // cout<<"2"<<phixgenome.substr(phixgenome.size()+pos-int(seq.size()),seq.size())<<endl;
		    //unsigned int posSeq=seq.size()-1;
		    pos++;//bring back to zero based
		    seqOut<<phixgenomeRC.substr(-1*pos,seq.size())<<endl;
		    // for(unsigned int posPhiX=(phixgenome.size()+pos-int(seq.size()));posPhiX<((phixgenome.size()+pos-int(seq.size()))+seq.size());posPhiX++,posSeq--){

		    // 	char b=complement(seq[posSeq]);
		    // 	//cout<<phixgenome[posPhiX]<<"\t"<<b<<endl;
		    // 	vecCov[posPhiX].cov++;
		    // 	if(phixgenome[posPhiX] != b &&
		    // 	   b                   != 'N'){
		    // 	    vecCov[posPhiX].base[ base2int( b) ]++;
		    // 	}
		    // }

		}else{
		    
		    //unsigned int posSeq=0;
		    pos--;//bring back to zero based
		    seqOut<<phixgenome.substr(pos,seq.size())<<endl;
		    // for(unsigned int posPhiX=pos;posPhiX<(pos+seq.size());posPhiX++,posSeq++){

		    // 	vecCov[posPhiX].cov++;
		    // 	if(phixgenome[posPhiX] != seq[posSeq] &&
		    // 	   seq[posSeq]         != 'N'){
		    // 	    vecCov[posPhiX].base[ base2int(seq[posSeq]) ]++;
		    // 	}
		    // }
		}
		// }
	    }
	    myFile.close();
	}else{
	    cerr << "Unable to open file "<<filename<<endl;
	    return 1;
	}
    }
    seqOut.close();
	
    return 0;
}

