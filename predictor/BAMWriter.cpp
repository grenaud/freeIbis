/*
 * BAMWriter
 * Date: Sep-04-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "BAMWriter.h"

BAMWriter::BAMWriter(string filename){
    myfile=new ofstream();
    
    myfile->open(filename.c_str(), ios::out | ios::binary);
    if (!myfile->is_open()) { 
	cerr<<"Cannot open file "<<filename<<endl;
	exit(1);
    }
}

BAMWriter::~BAMWriter(){
    close();
}

void BAMWriter::close(){
    myfile->close();    
}




inline unsigned char BAMWriter::bp2number(char b){
    //cerr<<"b "<<b<<endl;
    if(b == 'N')
	return valueN;
    if(b == 'A')
	return valueA;
    if(b == 'C')
	return valueC;
    if(b == 'G')
	return valueG;
    if(b == 'T')
	return valueT;
    cerr<<"Invalid bp "<<b<<endl;
    exit(1);
    
}




//void BAMWriter::writeZero(ofstream * myfile){
void BAMWriter::writeZero(){
    myfile->write((char*)&zero,sizeof(zero));
}


void BAMWriter::writeHeader(){ //(ofstream * myfile){
    //MAGIC FIELD
    const  unsigned char  stringHeader[] = "BAM\1";
    const unsigned int     lengthHeader  = 4;

    myfile->write((char*)&stringHeader,lengthHeader);
    
    //header + references, all zeros
    for(int i = 0;i<8; i++){
	//myfile->write((char*)&zero,sizeof(zero));
	writeZero();
    }
    
}




inline void BAMWriter::writeTag(//ofstream * myfile,
				string & tag){
    for(unsigned int i=0;i<tag.size();i++){
	unsigned char towrite=tag[i];
	myfile->write ((char*)&towrite,sizeof(towrite) );
    }
    writeZero();       
}



//subroutine to write an unmapped BAM record
inline void BAMWriter::writeSequence(//ofstream * myfile,
			  uint16_t flagToUse,
			  string & name,
			  string & firstSeq,
			  string & firstQual,
			  string & index1,
			  string & index1Q,
			  string & index2,
			  string & index2Q,
			  bool hasIndex1,
			  bool hasIndex2){

    unsigned int tagSize=0;

    if(hasIndex1 ){
	tagSize+= 2*(index1.size()+3+1); //3 for the XIZ and times two for the quality   +1 for the each zero byte (times 2: 1 for index, 1 for qual)

	if(index1.size() != index1Q.size()){
	    cerr<<"The size of the sequence for the first index and its associated quality differs"<<endl;
	    exit(1);
	}
    }

    if(hasIndex2){
	tagSize+= 2*(index2.size()+3+1); //3 for the XIZ and times two for the quality  +1 for the each zero byte (times 2:1 for index, 1 for qual)    

	if(index2.size() != index2Q.size()){
	    cerr<<"The size of the sequence for the second index and its associated quality differs"<<endl;
	    exit(1);
	}
    }

    if(firstSeq.size() != firstQual.size()){
	cerr<<"The size of the sequence and quality differs"<<endl;
	exit(1);
    }


    int coresize=8;


    unsigned int namelength=(name.size()+1);
    unsigned int queryLength=firstSeq.size();

    //COMPUTING BLOCK SIZE
    uint32_t sizeOfBlock= 
	coresize * sizeof(unsigned int)+  // core size
	namelength +                      //size of read name
	0 +                               // no cigar
	((queryLength+1) /2) +            //encoded sequence size
	firstQual.size() +                //size of qual
	(tagSize) ;                       //size of tag

    myfile->write ((char*)&sizeOfBlock, sizeof (sizeOfBlock));



    unsigned int maxunit=UINT_MAX;


    //8 fields of CORE
    //1 refID
    myfile->write ((char*)&maxunit,sizeof(maxunit));
    //2 pos
    myfile->write ((char*)&maxunit,sizeof(maxunit));
    //3 bin_mq_nl
    myfile->write ((char*)&namelength,sizeof(namelength));
    //4 flag_nc
    uint16_t flag=flagToUse;
    unsigned int flag_nc = flag<<16 | 0;
    myfile->write ((char*)&flag_nc,sizeof(flag_nc));    
    //5 l_seq
    myfile->write ((char*)&queryLength,sizeof(queryLength));    
    //6 next_refID
    myfile->write ((char*)&maxunit,sizeof(maxunit));
    //7 next_pos
    myfile->write ((char*)&maxunit,sizeof(maxunit));
    //8 t_len
    unsigned int zero_uint = 0;
    myfile->write ((char*)&zero_uint,sizeof(zero_uint) );

    //QUERY NAME
    for(unsigned int i=0;i<name.size();i++)
	myfile->write ((char*)&name[i],sizeof(unsigned char));
    unsigned char endofstring='\0';
    myfile->write ((char*)&endofstring,sizeof(unsigned char));



    //SEQ
    for(unsigned int i=0;i<(firstSeq.size()-1);i+=2){
	// cerr<<i<<endl;
	// cerr<<bp2number(firstSeq[i]) <<endl;
	// cerr<<bp2number(firstSeq[i+1])<<endl;
	unsigned char towrite=bp2number(firstSeq[i])<<4 | bp2number(firstSeq[i+1]);
	myfile->write ((char*)&towrite,sizeof(towrite) );
    }

    //final char in case of uneven seq
    if(firstSeq.size()%2 != 0){	
	unsigned char towrite=bp2number(firstSeq[ firstSeq.size()-1 ])<<4 | zero;
	myfile->write ((char*)&towrite,sizeof(towrite) );
    }


    //QUAL
    for(unsigned int i=0;i<firstQual.size();i++){
	unsigned char towrite=int(firstQual[i])-33;
	myfile->write ((char*)&towrite,sizeof(towrite) );
    }

    //tags
    string tag;
    if(hasIndex1 ){
	tag="XIZ"+index1;
	writeTag(tag);
	 tag="YIZ"+index1Q;
	writeTag(tag);
    }

    if(hasIndex2 ){
	tag="XJZ"+index2;
	writeTag(tag);
	tag="YJZ"+index2Q;
	writeTag(tag);
    }

 

}




void BAMWriter::writePairedSequence(//ofstream * myfile,
					   string name,
					   string firstSeq,
					   string firstQual,
					   string secondSeq,
					   string secondQual,
					   string index1,
					   string index1Q,
					   string index2,
					   string index2Q,
					   bool hasIndex1,
					   bool hasIndex2){
    writeSequence(//myfile,
		  flagFirstPair,
		  name,
		  firstSeq,
		  firstQual,
		  index1,
		  index1Q,
		  index2,
		  index2Q,
		  hasIndex1,
		  hasIndex2);

    writeSequence(//myfile,
		  flagSecondPair,
		  name,
		  secondSeq,
		  secondQual,
		  index1,
		  index1Q,
		  index2,
		  index2Q,
		  hasIndex1,
		  hasIndex2);

}

void BAMWriter::writeSingleSequence(//ofstream * myfile,
				     string name,
				     string firstSeq,
				     string firstQual,
				     string index1,
				     string index1Q,
				     string index2,
				     string index2Q,
				     bool hasIndex1,
				     bool hasIndex2){
    writeSequence(//myfile,
		flagSingleReads,
		  name,
		  firstSeq,
		  firstQual,
		  index1,
		  index1Q,
		  index2,
		  index2Q,
		  hasIndex1,
		  hasIndex2);



}
