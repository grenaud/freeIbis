#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <map>
#include <vector>
#include <math.h>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h> 
#include <sys/time.h>

#include "func.h"

using namespace std;

int maxEditDist         = 5;
int ctrlSeqINDEXlength  = 7;
int seedlength          = 30;
int wordlength          = 10;
int qualOffset          = 33;

typedef struct { 
    char baseC;
    int  qualC;
} base;


inline int hammingDistance(const string & s1,
			   const string & s2){

    if(s1.size() != s2.size()){ 	
	cerr<<"ERROR: the hammingDistance() function cannot be called for strings of different lengths"<<endl; 	
	exit(1);     
    }

    int substitutions=0;
    for(unsigned int i=0;i<s1.size();i++){    
	// if(s1[i] == 'N' ||
	//    s2[i] == 'N' )
	//     continue;
	   

	if(s1[i] != s2[i]){
	    substitutions++;
	}
    }

    return substitutions;
}

inline bool isWithinHammingDistance(const string & s1,
			   const string & s2,
			   int maxDist){
    if(s1.size() != s2.size()){ 	
	cerr<<"ERROR: the isWithinHammingDistance() function cannot be called for strings of different lengths"<<endl; 	
	exit(1);     
    }

    int substitutions=0;
    for(unsigned int i=0;i<s1.size();i++){    
	if(s1[i] != s2[i]){
	    substitutions++;
	}
	if(substitutions > maxDist){
	    return false;
	}
    }
    return true;
}



inline base bin2base(unsigned char c ){
    base toreturn;
    toreturn.baseC="ACGT"[c % 4];
    c /= 4;
    toreturn.qualC=(c % 64);
    if(toreturn.qualC == 33){
    	toreturn.baseC='N';
    }
    return toreturn;
}


inline bool isFile(const string & dir){
    struct stat st;
    if(stat(dir.c_str(),&st) == 0)
	if( st.st_mode & S_IFREG )
	    return true;
    
    return false;
}


inline unsigned int hashword(const string & dnaString){
    unsigned int toReturn = 0;

    for(int i=0;i<( dnaString.size() );i++){
	toReturn = toReturn<<2;

	if(dnaString[i] == 'A'){
	    //nothing to do
	    continue;

	}

	if(dnaString[i] == 'C'){
	    toReturn |= 1;
	    continue;
	}

	if(dnaString[i] == 'G'){
	    toReturn |= 2;
	    continue;
	}

	if(dnaString[i] == 'T'){
	    toReturn |= 3;
	    continue;
	}   
 
	if(dnaString[i] == 'N'){
	    toReturn |= rand()%4;
	    continue;	    
	}

	cerr<<"hashword(): invalid nucleotide in string: "<<dnaString<<endl;
	exit(1);    

    }
    return toReturn;
}

template <typename T>
string stringify(const T i){
    stringstream s;
    s << i;
    return s.str();
}

int main (int argc, char *argv[]) {
    timeval time;
    gettimeofday(&time, NULL);
    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );

    const string usage=string(argv[0])+"\t [Directory with basecalls] [PhiX fasta file]"+
	"Options:\n"+
	"\tMandatory\n"+
	"\t\t-f [First cycle PhiX]\n"+
	"\t\t-n [Number of cycle for PhiX]\n"+
	"\t\t-l [lane#]\n"+
	"\t\t-t [tile#]\n"+
	"\tOptional:\n"+
	"\t\t-i [First cycle index]\n"+
	"\t\t-s [Index sequence(s)]\tComma separated if multiple indices were used\n"+
	"\t\t-d [distance (Default "+stringify( maxEditDist )+")]\n"+
	"\t\t-w [word length (Default "+stringify( wordlength )+")]\n\n"+

	"Example:\n" +
	"\tbcl2phix -f 1  -i 77 -l 1 -s CGATTCG -t 1101 /mnt/solexa/120217_SOLEXA-GA02_00048_PEdi_AW_CB_MM/Data/Intensities/BaseCalls/ /mnt/solexa/Genomes/phiX/control/whole_genome.fa\n";
    //CGATTCG

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<usage<<endl;
	
	return 1;
    }

    string lane             = "";
    string tile             = "";
    int firstCycle          = -1;
    int numberCycle         = -1;
    int lastCycle           = -1;

    int firstCycleIDX       = -1;
    string idxSequenceTemp           = "";
    vector<string> idxSequence;

    bool   idxSequenceSpecified  = false;

    bool noIndex= true;    
    for(int i=1;i<(argc-1);i++){ //all but the last arg

	if(string(argv[i]) == "-f" ){
	    firstCycle=atoi(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-n" ){
	    numberCycle=atoi(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-d" ){
	    maxEditDist=atoi(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-w" ){
	    wordlength=atoi(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-i" ){
	    firstCycleIDX=atoi(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-l" ){
	    lane=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-t" ){
	    tile=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-s" ){
	    idxSequenceTemp=string(argv[i+1]);
	    idxSequenceSpecified=true;
	    i++;
	    continue;
	}


    }

    if(idxSequenceSpecified){
	size_t lastfound=-1;
	while(true){
	    size_t found = idxSequenceTemp.find(',',lastfound+1);
	    if (found!=string::npos){
		idxSequence.push_back(idxSequenceTemp.substr(lastfound+1,found-lastfound-1));
		lastfound=found;
	    }else{
		idxSequence.push_back(idxSequenceTemp.substr(lastfound+1));
		break;
	    }
	}

	ctrlSeqINDEXlength = int(idxSequence[0].size());
	for(unsigned int i=1;i<idxSequence.size();i++){
	    if(int(idxSequence[i].size()) != ctrlSeqINDEXlength){
		cerr<<"Index= "<<idxSequence[i]<<" does not have the same length as the other "<<endl;
		return 1;
	    }
	}
	           
    }

    if(lane.empty()){
	cerr<<"lane must be defined"<<endl;
	return 1;
    }

    if(tile.empty()){
	cerr<<"tile must be defined"<<endl;
	return 1;
    }

    if(firstCycle == -1){
	cerr<<"first cycle of must be defined"<<endl;
	return 1;
    }

    if(numberCycle == -1){
	cerr<<"number of cycles of must be defined"<<endl;
	return 1;
    }

    if(firstCycleIDX != -1){
	noIndex=false;
	if(idxSequence.empty()){
	    cerr<<"sequence for the index must be defined"<<endl;
	    return 1;
	}

    }

    string directory        = string(argv[argc-2]);
    string phixFileName     = string(argv[argc-1]);
    lastCycle = firstCycle + numberCycle;




    //reading fasta file 
    string line;
    ifstream phixFile;
    string phixgenome   = "";
    string phixgenomeRC = "";
    int * hashWord2Pos;
    try{
	hashWord2Pos  = new int[ int(pow(2,2*wordlength))];
    }catch(bad_alloc&) {
	cerr<<"Failed to allocate memory for a hash of  "<<int(pow(2,2*wordlength))<<" elements"<<endl;
	return 1;
    }

    for(int i=0;i<int(pow(2,2*wordlength));i++){
	hashWord2Pos[i] = 0;
    }



    phixFile.open(phixFileName.c_str(), ios::in);
    
    if (phixFile.is_open()){
	getline (phixFile,line);
	while ( getline (phixFile,line)){
	    transform(line.begin(), line.end(), line.begin(), ::toupper);
	    phixgenome+=line;
	}
	phixFile.close();
    }else{
	cerr << "Unable to open file "<<phixFileName<<endl;
	return 1;
    }


    phixgenomeRC=reverseComplement(phixgenome);

    map<string,unsigned int> word2count;
    map<string,int> word2pos; //1 based (pos=+ strand,neg=- strand)

    for(int i=0;i<=(phixgenome.size()-wordlength);i++){
	// cout<<phixgenome.substr(i,wordlength)<<endl;
	
	word2pos[ phixgenome.substr(i,wordlength) ]  = i+1;
	
	if(word2count.find(phixgenome.substr(i,wordlength)) == word2count.end()){
	    word2count[ phixgenome.substr(i,wordlength) ] = 1;
	}else{
	    word2count[ phixgenome.substr(i,wordlength) ] ++;
	}
    }

    for(int i=0;i<=(phixgenomeRC.size()-wordlength);i++){
    	//cout<<phixgenomeRC.substr(i,wordlength)<<endl;
	word2pos[ phixgenomeRC.substr(i,wordlength) ]  = -(i+1);
    	if(word2count.find(phixgenomeRC.substr(i,wordlength)) == word2count.end()){
    	    word2count[ phixgenomeRC.substr(i,wordlength) ] = 1;
    	}else{
    	    word2count[ phixgenomeRC.substr(i,wordlength) ] ++;
    	}	


    }
    

    for(map<string,unsigned int>::iterator iter = word2count.begin(); 
	iter != word2count.end(); iter++) {
	// cout<<iter->first<<"\t"<<iter->second<<"\t"<<word2pos[ iter->first ] <<endl;

	if(iter->second == 1){//unique hit
	    //cerr<<"uniq "<<iter->first<<endl;
	    hashWord2Pos[ hashword(iter->first) ] = word2pos[ iter->first ];
	}else{
	    //cerr<<"not uniq "<<iter->first<<"\t"<<iter->second<<endl;
	}
    }
    // cout<<phixgenome<<endl;
    // return 1;












    ///////////////////////////////////////
    //      BEGIN READING BCL FILES      //
    ///////////////////////////////////////

    stringstream ss;
    fstream mybclfileIndex [firstCycleIDX+ctrlSeqINDEXlength-1];
    fstream mybclfilePhix  [numberCycle];

    unsigned int numberOfClustersFirst=0 ;
    // unsigned int numberOfCtrlClusters =0 ;

    
    //open index files
    if(!noIndex){
	for(int cycle=firstCycleIDX;cycle<=(firstCycleIDX+ctrlSeqINDEXlength-1);cycle++){
	    ss.str(std::string());
	    ss<<cycle;
	    string bclFile=directory+"/L00"+lane+"/C"+ ss.str()+".1/s_"+lane+"_"+tile+".bcl";
	    //cerr<<bclFile<<endl;
	    if(!isFile(bclFile)){
		cerr<<"Unable to find file "<<bclFile<<endl;
		return 1;
	    }
	    mybclfileIndex[cycle-firstCycleIDX].open(bclFile.c_str(),ios::in|ios::binary);
	    if (!mybclfileIndex[cycle-firstCycleIDX]) {
		cerr<<"Unable to read index file "<<bclFile<<endl;
		return 1;
	    }
	    unsigned int numberOfClusters ;
	    mybclfileIndex[cycle-firstCycleIDX].read((char*)&numberOfClusters, sizeof (int));
	    //cerr<<"numberOfClusters "<<numberOfClusters<<endl;
	    if(numberOfClustersFirst == 0 ){
		numberOfClustersFirst=numberOfClusters;
	    }else{
		if(numberOfClustersFirst!=numberOfClusters){
		    cerr<<"Number of clusters in "<<bclFile<<" is different from the other cycles, exiting\n";
		    return 1;
		}
	    }
	}
    }

    for(int cycle=firstCycle;cycle<=(firstCycle+numberCycle-1);cycle++){
	ss.str(std::string());
	ss<<cycle;
	string bclFile=directory+"/L00"+lane+"/C"+ ss.str()+".1/s_"+lane+"_"+tile+".bcl";
	//cerr<<bclFile<<endl;
	if(!isFile(bclFile)){
	    cerr<<"Unable to find file "<<bclFile<<endl;
	    return 1;
	}


	mybclfilePhix[cycle-firstCycle].open(bclFile.c_str(),ios::in|ios::binary);

	if (!mybclfilePhix[cycle-firstCycle]) {
	    cerr<<"Unable to read file "<<bclFile<<endl;
	    return 1;
	}

	unsigned int numberOfClusters ;
	mybclfilePhix[cycle-firstCycle].read((char*)&numberOfClusters, sizeof (int));
	//cerr<<"numberOfClusters "<<numberOfClusters<<endl;
	if(numberOfClustersFirst == 0 ){
	    numberOfClustersFirst=numberOfClusters;
	}else{
	    if(numberOfClustersFirst!=numberOfClusters){
		cerr<<"Number of clusters in "<<bclFile<<" is different from the other cycles, exiting\n";
		return 1;
	    }
	}
    }

    int firstCycleToPrint  = firstCycle-1;
    int lastCycleToPrint   = lastCycle-1;

    for(unsigned int clusterN=0;clusterN<numberOfClustersFirst;clusterN++){
	//cout<<"cluster #"<<i<<endl;
	if(!noIndex){
	    string sequenceIDX="";
	    string qualityIDX ="";

	    for(int idxcycle=0;idxcycle<ctrlSeqINDEXlength;idxcycle++){
		char toread;
		mybclfileIndex[idxcycle].read(&toread, sizeof (char));
		base returned=bin2base(toread);
		sequenceIDX   +=        returned.baseC;
		qualityIDX    +=char(33+returned.qualC);
	    }

	    //if(sequenceIDX != idxSequence){
	    //	for(unsigned int i=0;i<int(idxSequence.size());i++){
	    bool oneMatch=false;
	    for(unsigned int idxIdxSeq=0;idxIdxSeq<idxSequence.size();idxIdxSeq++){
		if(isWithinHammingDistance(sequenceIDX,idxSequence[idxIdxSeq],1)){ // is within hamming dist, set to true and break
		    oneMatch=true;
		    break;
		}
	    }

	    if(!oneMatch){ // if not a single one matched

		//advance all fps
		char toread;
		for(int phixcycle=0;phixcycle<numberCycle;phixcycle++){
		    mybclfilePhix[phixcycle].read(&toread, sizeof (char));
		}
		continue;
	    }
	    // cout<<sequenceIDX<<endl;
	}


	string sequence  = "";
	string quality   = "";

	int maxQual      =  0;
	int maxQualIndex =  0;
	int scdMaxQual      =  0;
	int scdMaxQualIndex =  0;


	int  currentQual =  0;
	int lastQuals[wordlength];

	//First cycle
	int phixcycle=0;
	char toread;
	mybclfilePhix[phixcycle].read(&toread, sizeof (char));
	
	base returned=bin2base(toread);
	sequence   +=returned.baseC;
	quality    +=char(33+returned.qualC);
	if(phixcycle<wordlength){
	    lastQuals[0] = returned.qualC;
	}
	//remaining cycles
	for(int phixcycle=1;phixcycle<numberCycle;phixcycle++){
	    
	    mybclfilePhix[phixcycle].read(&toread, sizeof (char));

	    returned    = bin2base(toread);
	    sequence   += returned.baseC;
	    quality    += char(33+returned.qualC);

	    if(phixcycle<wordlength){//before wordlength
		//lastFirstQual = 
		currentQual          += returned.qualC;
		lastQuals[phixcycle]  = returned.qualC;
		maxQual               = currentQual;
	    }else{
		if(maxQual < currentQual){ //new qual
		    scdMaxQual      = maxQual;
		    scdMaxQualIndex = maxQualIndex;

		    maxQual      = currentQual;
		    maxQualIndex = phixcycle-wordlength;
		}else{
		    if(scdMaxQual < currentQual){ //new second max qual
			scdMaxQual      = currentQual;
			scdMaxQualIndex = phixcycle-wordlength;
		    }
		}

		currentQual                     -= lastQuals[phixcycle%wordlength]; //returned.qualC;
		currentQual                     += returned.qualC;
		lastQuals[phixcycle%wordlength]  = returned.qualC;
	    }
	    // cout<<phixcycle<<"\t"<<currentQual<<endl;
	}
	
	// for(unsigned int k=0;k<quality.size();k++){
	//     cout<<(int(quality[k])-qualOffset)<<endl;
	// }
	int indexPhixFound;

	unsigned int indexInHash = hashword( sequence.substr(maxQualIndex,wordlength) );


	// unsigned int indexInHash = 0; //hashword( sequence.substr(maxQualIndex,wordlength) );

	// for(int i=maxQualIndex;i<( maxQualIndex+wordlength );i++){
	//     indexInHash = indexInHash<<2;
	//     //cout<<sequence[i]<<endl;
	//     if(sequence[i] == 'A'){
	// 	//nothing to do
	// 	continue;		
	//     }
	    
	//     if(sequence[i] == 'C'){
	// 	indexInHash |= 1;
	// 	continue;
	//     }
	    
	//     if(sequence[i] == 'G'){
	// 	indexInHash |= 2;
	// 	continue;
	//     }
	    
	//     if(sequence[i] == 'T'){
	// 	indexInHash |= 3;
	// 	continue;
	//     }   
	    
	//     //invalid nucleotide in string
	//     goto endloop;	    
	// }

	indexPhixFound = hashWord2Pos[indexInHash];
	// cout<<"reached"<<indexPhixFound<<"\t"<<indexInHash<<"\t"<<hashword( sequence.substr(maxQualIndex,wordlength)) <<endl;
	// if(clusterN == 3284){
	//     cout<<sequence<<endl;
	//     cout<<quality<<endl;
	//     cout<<sequence.substr(maxQualIndex,wordlength)<<endl;

	//     // cout<<phixgenomeRC.substr( st ,numberCycle)<<endl;
	//     cout<<indexPhixFound<<endl;
	//     int st = (indexPhixFound-maxQualIndex-1);
	//     cout<<phixgenome.substr( st ,numberCycle)<<endl;
	//     cout<<hammingDistance(sequence,phixgenome.substr( st ,numberCycle))<<endl;
	//     cout<<endl;
	//     return 1;
	// }

	if(indexPhixFound != 0){
	    int st;
		
	    if(indexPhixFound<0){//- strand		
		st = (-1*indexPhixFound-maxQualIndex-1);
		if(st<0)
		    continue;
		if( (st+numberCycle) >= phixgenome.size())
		    continue;
		
		if(isWithinHammingDistance(sequence,phixgenomeRC.substr( st ,numberCycle),maxEditDist)){
		    cout<<lane<<"\t"<<tile<<"\t"<<clusterN<<"\tIDX\t"<<firstCycleToPrint<<"\t"<<lastCycleToPrint<<"\t"<<sequence<<"\t"<<( -1*(st+1) )<<endl;
		}

	    }else{               //+ strand 		
		st = (indexPhixFound-maxQualIndex-1);
		if(st<0)
		    continue;
		if( (st+numberCycle) >= phixgenome.size())
		    continue;
		
		if(isWithinHammingDistance(sequence,phixgenome.substr( st ,numberCycle),maxEditDist)){
		    cout<<lane<<"\t"<<tile<<"\t"<<clusterN<<"\tIDX\t"<<firstCycleToPrint<<"\t"<<lastCycleToPrint<<"\t"<<sequence<<"\t"<<(st+1)<<endl;
		}

	    }
	    
	}else{
	    
	    if(	 scdMaxQualIndex != maxQualIndex){//try second highest 

		indexInHash = hashword( sequence.substr(scdMaxQualIndex,wordlength) );
		indexPhixFound = hashWord2Pos[indexInHash];

		if(indexPhixFound != 0){//found
		    int st;
		
		    if(indexPhixFound<0){//- strand		
			st = (-1*indexPhixFound-scdMaxQualIndex-1);
			if(st<0)
			    continue;
			if( (st+numberCycle) >= phixgenome.size())
			    continue;
		
			if(isWithinHammingDistance(sequence,phixgenomeRC.substr( st ,numberCycle),maxEditDist)){
			    cout<<lane<<"\t"<<tile<<"\t"<<clusterN<<"\tIDX\t"<<firstCycleToPrint<<"\t"<<lastCycleToPrint<<"\t"<<sequence<<"\t"<<( -1*(st+1) )<<endl;
			}

		    }else{               //+ strand 		
			st = (indexPhixFound-scdMaxQualIndex-1);
			if(st<0)
			    continue;
			if( (st+numberCycle) >= phixgenome.size())
			    continue;
		
			if(isWithinHammingDistance(sequence,phixgenome.substr( st ,numberCycle),maxEditDist)){
			    cout<<lane<<"\t"<<tile<<"\t"<<clusterN<<"\tIDX\t"<<firstCycleToPrint<<"\t"<<lastCycleToPrint<<"\t"<<sequence<<"\t"<<(st+1)<<endl;
			}

		    }
	    
		}else{
		    //give up
		    // cerr<<"notfound\t"<<clusterN<<"\t"<<sequence<<"\t"<<sequence.substr(maxQualIndex,wordlength)<<"\t"<<sequence.substr(scdMaxQualIndex,wordlength)<<endl;
		}
	    }


	}
	///cout<<endl<<endl;
	//return 1;
	// if(sequence == ctrlSeqINDEX){
	//     numberOfCtrlClusters++;
	// }
	// endloop:
	// 	continue;	
    }//end for over each clusters


    ///////////////////////////////////////
    //      END READING BCL FILES      //
    ///////////////////////////////////////

    // cout<<lane<<"\t"<<tile<<"\t"<<numberOfCtrlClusters<<"\t"<<numberOfClustersFirst<<"\t"<<100.0*(double(numberOfCtrlClusters)/double(numberOfClustersFirst))<<"%"<<endl;

    return 0;
}

