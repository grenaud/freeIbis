// :vim:ts=8:
//
// This wants to be streamlined.  Generally speaking, we want to
// minimize disk I/O, especially excessive seeking.  Depending on the
// layout of the input data, different strategies make sense.  We'd like
// to handle these formats for position information:
//
// - Firecrest format (text, one position per line, complete with
//   intensities)
// - IPAR format (text, one position per line)
// - MiSeq format (locs file, one position per record)
// - HiSeq format (clocs file, one position per record)
//
// and these intensity formats:
//
// - Firecrest (text, intensities in columns)
// - IPAR (text, one block per cycle, one line per cluster with four
//   intensities)
// - CIF (binary, one file per cycle, one block per channel)
//
// Output is BAM, no discussion about it.
//
// For Firecrest, we can just loop over the input, go record by record
// and produce BAM directly.  Both the calling loop over records and the
// compression loop over output chunks could be parallelized, but we
// don't really care.
//
// For CIF, we must go cycle by cycle, in big blocks to cut down on
// seeking.  The output needs to be buffered before creating BAM.  (We
// could buffer input instead, but it's much larger.) Suppose we're on
// HiSeq, we get 24MB of intensities per cycle, which will result in 3MB
// of basecalls, or 750MB of base calls for the whole tile over all
// cycles.  I think this is still practical, and remains practical for
// anticipated read lengths and cluster densities.  Unless Illumina
// dumps way bigger files on us, we're not going to run out of memory.
// However, we need to parallelize internally; calling multiple
// instances of the predictor isn't going to cut it.
//
// So, the strategy goes like this:  A small number of input threads
// fires off requests to read a few CIF files ahead.  A larger number of
// worker threads calls small blocks of bases for the current cycle.
// When done with some cycle, that CIF file can be freed and we continue
// on the next cycle.  When we're done with the last cycle, we have an
// enormous block of base calls (base + quality).  We can roughly
// predict the size of the resulting BAM, so we chop this hunk into
// blocks, and a number of worker threads can convert a block at a time
// to BAM, gzip it, and send it to output.  Effectively, we have two
// very separate steps: reading and calling first, then compressing and
// writing BAM.  It doesn't look as if overlapping those two processes
// for different tiles is going to be a big win.
//
// For IPAR, we read the intensities into memory.  This was never a big
// problem, but is inelegant nonetheless.  We could instead treat IPAR
// files like CIF, read cycle by cycle and stream them.  Might be
// worthwhile, especially if it cleans the code up.
//
//
// Right now, we're not going to read CIF files in chunks.  It's
// complicated, causes more seeking, and is not necessary for the
// foreseeable future.  If those files get much bigger, we can split
// them.  The same if even more true for IPAR files.
//
// Reading from locs files could be parallelized, if we wanted to.
// Unfortunately, this is not the case for the more interesting clocs
// files; those are constructed in such a way that they must be read
// sequentially.  So we will necessarily have only one thread reading the
// locations.  This should still be fine, creating and compressing BAM
// should be more computational work, and the latter can be
// parallelized.


extern "C" {
#include "../libocas_v093/lib_svmlight_format.h"
#include "../libocas_v093/libocas.h"
}

#include "input.h"
#include "erf.h"
#include "BAMWriter.h"

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>
#include <libgen.h>
#include <stdint.h>
#include <sys/stat.h>

//Max length of a modelfile line
#define MODELFILE_MAXLINELEN 1000000

//for strerror, we should not have these 3:
#include <errno.h>
#include <stdio.h>
#include <string.h>

using namespace std;

char qscore_offset = 33;
double max_qscore = 40;
double min_qscore = 1;

static int cyclesInModelFile ; //to count the # of cycles in the model file
static unsigned int numberOfClustersFound; //# of cluster found in intensity files
int outformat = 0;
int coordtype = 3;
bool recalibrateQuals;
bool quick;
BAMWriter * bw;

typedef struct {
    double k;
    double interceptPreK;
    double interceptPostK;
    double slopePreK;
    double slopePostK;
}  LinearRegParams;

typedef struct {
    double * wMatrix;
    int nDimensions; //nDim
    int nColumns; //  nY or nCols;
}  MODELDATA;


typedef struct {
    //    int label;
    int pred_label;
    double * scores;
}  RETURNEDRESULTS;

struct Tmodel {
    // STRUCTMODEL model;
    // STRUCT_LEARN_PARM sparm;
    MODELDATA * matrix;

    MODELDATA       ** qualScoreSVM;
    LinearRegParams ** linearRegParams;

    char type;
    
    double factors[4][4];
    double sfactors[4][4];
    double rclassifier[4][4];

    bool use_factor;
};


// void free_struct_model(Tmodel& m
struct FreeModel { 
    void operator()( Tmodel& m ) { 
	free(m.matrix->wMatrix);	    	
	free(m.matrix);	    

	if(recalibrateQuals){    
	    for(int qssvm=0;qssvm<4;qssvm++){		
		free(m.qualScoreSVM[ qssvm ]->wMatrix);
		delete(m.linearRegParams[ qssvm ]);
		free(m.qualScoreSVM[ qssvm ]);
	    }
	}

	delete[] m.linearRegParams ;
	delete[] m.qualScoreSVM ;
	
    } ;
} ;

string intStringify(int i) {
    stringstream s;
    s << i;
    return s.str();
}

// inline void deleteRETURNEDRESULTS(RETURNEDRESULTS * todelete){
//     delete(todelete->scores);
//     delete(todelete);
// }

inline static unsigned char char2bin(char b){
    if(b == 'A')
	return 0;
    if(b == 'C')
	return 1;
    if(b == 'G')
	return 2;
    if(b == 'T')
	return 3;
    return 0;
}

#if 0
double density(double x,double mu, double sigma)
{
    double y = (x-mu)/sigma;
    return (1.0/(sigma*sqrt2pi))*exp((y*y)/-2.0);
}
#endif

void append_zero_quality(string &s, string &qb, string &qa, string &qc, string &qg, string &qt )
{
    s.push_back( 'N' ) ;
    qb.push_back( qscore_offset ) ;
    qa.push_back( qscore_offset ) ;
    qc.push_back( qscore_offset ) ;
    qg.push_back( qscore_offset ) ;
    qt.push_back( qscore_offset ) ;
}


void append_seq_quality(string &s, string &qb, string &qa, string &qc, string &qg, string &qt, const RETURNEDRESULTS * y, const Tmodel* model,ofstream * toWrite=NULL)
{
    int call_class = y->pred_label-1; // CLASS DETERMINED BY SVM LIGHT

    double qualities[5] = {0,0,0,0,0};

    if(recalibrateQuals){
	double bestDist =-1000.0;
	int calledDist  =-1;
	for (unsigned int cclass = 0; cclass < 4; cclass++){

	    MODELDATA * matrixToUse;
	    if(!quick)
		matrixToUse=model->qualScoreSVM[cclass];

	    double sumForClass=0.0;
	    double sumForAll=0.0;
	    if(quick){
		for( unsigned oclass = 0 ; oclass < 4 ; oclass++ ){
		    if(cclass != oclass){
			sumForClass+=(1.0/(1.0+exp(-1.0*y->scores[oclass])));
		    }
		    sumForAll+=(1.0/(1.0+exp(-1.0*y->scores[oclass])));
		}
		double postprob=log( sumForClass / sumForAll );
		if(postprob < model->linearRegParams[cclass]->k){
		    qualities[cclass]= postprob*model->linearRegParams[cclass]->slopePreK  + model->linearRegParams[cclass]->interceptPreK;
		}else{
		    qualities[cclass]= postprob*model->linearRegParams[cclass]->slopePostK + model->linearRegParams[cclass]->interceptPostK;
		}
	    }else{
		for( unsigned oclass = 0 ; oclass < 4 ; oclass++ ){
		    sumForClass+=y->scores[oclass]*matrixToUse->wMatrix[LIBOCAS_INDEX(oclass,0,matrixToUse->nDimensions)];     
		}
		sumForClass+=matrixToUse->wMatrix[LIBOCAS_INDEX(4,0,matrixToUse->nDimensions)];     

		if(sumForClass < model->linearRegParams[cclass]->k){
		    qualities[cclass]= sumForClass*model->linearRegParams[cclass]->slopePreK  + model->linearRegParams[cclass]->interceptPreK;
		}else{
		    qualities[cclass]= sumForClass*model->linearRegParams[cclass]->slopePostK + model->linearRegParams[cclass]->interceptPostK;
		}
	    }


	    if(y->scores[cclass]  > bestDist){
		bestDist=y->scores[cclass];
		calledDist=cclass;
	    }
	}

	if(call_class >= 0)
	    calledDist=call_class;
	// cout<<"calledDist "<<calledDist<<endl;

	MODELDATA * matrixToUse;;
	if(!quick)
	    matrixToUse=model->qualScoreSVM[calledDist];
	double sumForClass=0.0;
	double sumForAll=0.0;

	if(quick){
	    for( unsigned oclass = 0 ; oclass < 4 ; oclass++ ){
		if(call_class != int(oclass)){
		    sumForClass+=(1.0/(1.0+exp(-1.0*y->scores[oclass])));//y->scores[oclass];
		}
		sumForAll+=(1.0/(1.0+exp(-1.0*y->scores[oclass])));
	    }
	    double postprob=log( sumForClass / sumForAll );

	    if(postprob < model->linearRegParams[calledDist]->k){
		qualities[4]= postprob*model->linearRegParams[calledDist]->slopePreK  + model->linearRegParams[calledDist]->interceptPreK;
	    }else{
		qualities[4]= postprob*model->linearRegParams[calledDist]->slopePostK + model->linearRegParams[calledDist]->interceptPostK;
	    }
	}else{
	    for( unsigned oclass = 0 ; oclass < 4 ; oclass++ ){
		sumForClass+=y->scores[oclass]*matrixToUse->wMatrix[LIBOCAS_INDEX(oclass,0,matrixToUse->nDimensions)];     
	    }
	    sumForClass+=matrixToUse->wMatrix[LIBOCAS_INDEX(4,0,matrixToUse->nDimensions)];     
	    if(sumForClass < model->linearRegParams[calledDist]->k){
		qualities[4]= sumForClass*model->linearRegParams[calledDist]->slopePreK  + model->linearRegParams[calledDist]->interceptPreK;
	    }else{
		qualities[4]= sumForClass*model->linearRegParams[calledDist]->slopePostK + model->linearRegParams[calledDist]->interceptPostK;
	    }
	}

    }else{
	//    if ((model->use_factor) && (outformat > 0) && (call_class != -1)){
	if ( (outformat == 1 ) && (call_class != -1) ){ //4q
	    double numerator[4];
	    double denominator = log(0.0) ;
	    for (unsigned int cclass = 0; cclass < 4; cclass++)	{
		numerator[cclass] = log(0.0) ;
		for( unsigned oclass = 0 ; oclass != 4 ; ++oclass )
		    add_logs( numerator[cclass],
			      log( model->rclassifier[cclass][oclass] ) +
			      log_cum_density(
					      //y.scores[oclass+1], under svmlight
					      y->scores[oclass],
					      model->factors[cclass][oclass],
					      model->sfactors[cclass][oclass] ) ) ;

		add_logs( denominator, numerator[cclass] ) ;
	    }

	    for (unsigned int cclass = 0; cclass < 4; cclass++)
		qualities[cclass] = -10.0/log(10.0) * (numerator[cclass]-log(3)-denominator); 

	    double qmin = qualities[0];
	    call_class = 0 ;
	    for (unsigned int cclass = 1; cclass < 4; cclass++)	{
		if (qualities[cclass] < qmin) { qmin = qualities[cclass]; call_class = cclass; }
	    }

	    qualities[4] = qmin ;
	}

	//if ((model->use_factor) && (call_class != -1)){
	if ( (call_class != -1) ){
	    double odensity = log(0.0), density = log(0.0) ;
	    for (int cclass = 0; cclass < 4; cclass++){
		double z =
		    log( model->rclassifier[call_class][cclass] ) +
		    log_cum_density(
				    // y.scores[cclass+1], under svmlight
				    y->scores[cclass],
				    model->factors[call_class][cclass],
				    model->sfactors[call_class][cclass] ) ;

		add_logs( density, z ) ;
		if (cclass != call_class) add_logs( odensity, z ) ;
	    }
	    qualities[4] = -10.0/log(10.0) * ( odensity - density ) ;
	}

    }


    string* qx[] = { &qa, &qc, &qg, &qt, &qb } ;
    bool good = false ;
    if(quick)
	good = true ;
    else
	for( size_t i = 0 ; i != 4 ; ++i )
	    good |= y->scores[i] >= 0 ;

    if( good ){
        for( size_t i = 0 ; i != 5 ; ++i ){
            qx[i]->push_back( char(qscore_offset + round( max(min_qscore,min( max_qscore,qualities[i])) ) ));
        }
        s.push_back( "NACGT"[1+call_class] ) ;
    }else{
        for( size_t i = 0 ; i != 5 ; ++i )
            qx[i]->push_back( qscore_offset ) ;
        s.push_back( 'N' ) ;
    }
};

void printModel(const MODELDATA * model){
    cout<<"Number of dims "<<model->nDimensions<<endl;
    cout<<"Number of cols "<<model->nColumns<<endl;

    for(int j=0; j < model->nColumns; j++) {
	for(int i=0; i < model->nDimensions; i++) {
	    cout<<"model ["<<i<<"]["<<j<<"]="<<model->wMatrix[LIBOCAS_INDEX(i,j,model->nDimensions)]<<endl;
	}
    }

}


LinearRegParams * readLinearParams(string inputParamsFile){
    string line;
    LinearRegParams * toReturn =new LinearRegParams;
    int lineCounter=0;
    ifstream valuesin;

    valuesin.open(inputParamsFile.c_str(), ios::in);  // open the streams  

    if (valuesin.is_open()) {  
	while ( getline (valuesin,line) ){   
	    istringstream stringStr( line ) ;    
	    double tempor;
	    if( stringStr >> tempor ) { 
	    }else{    
		cerr<<"Verify the value of k "<<inputParamsFile<<" wrong line #"<<lineCounter<<" line=#"<<line<<"#"<<endl;    
		exit(1);    
	    }	

	    if(lineCounter==0)
		toReturn->k=tempor;
	    if(lineCounter==1)
		toReturn->interceptPreK=tempor;
	    if(lineCounter==2)
		toReturn->slopePreK=tempor;
	    if(lineCounter==3)
		toReturn->interceptPostK=tempor;
	    if(lineCounter==4)
		toReturn->slopePostK=tempor;


	    lineCounter++;
	}     
	valuesin.close();    
    }else{
	cout << "Unable to open file="<<inputParamsFile<<endl;
	exit(1);
    }

    return toReturn;
}


/** 

   This function reads a model file created by 
   msvmocas and returns the address of a MODELDATA 
   struct created on the heap.

   This looks like C code because it was taken from linclass.c
   in the libocas package.
	 
   @param[in]     _inputMatrixFile Pointer to an array of char with the filename
   @return Pointer to a MODELDATA object created on the heap

*/
MODELDATA * readW(const char * inputMatrixFile){
    FILE *fid;         //file pointer
    char *line;        //temp var for chars
    uint32_t nLines = 0;    // number of lines
    char *endptr, *begptr;
    uint32_t i, j;
    double val;
    uint32_t nCols = 0, tmp_nCols;
    int go = 1;
    double * wLocal;
    MODELDATA * toReturn;

    toReturn = (MODELDATA *)calloc(1, sizeof(MODELDATA));

    line = (char *)calloc(MODELFILE_MAXLINELEN, sizeof(char));
    if( line == NULL ) {
	fprintf(stderr,"Not enough memmory to allocate line buffer.\n");
	exit(1);
    }

    fid = fopen(inputMatrixFile, "r");
    if(fid == NULL) {
	fprintf(stderr,"Cannot open model file. %s\n",inputMatrixFile);
	exit(1);
    }
  
    /* read the first line */
    if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ){
	fprintf(stderr,"Empty example file for %s\n",inputMatrixFile);
	exit(1);
    }else{
	nLines = 1;
	begptr = line;
	while(1){
	    val = strtod(begptr, &endptr);

	    if(val == 0 && begptr == endptr)
		break;

	    nCols++;
	    begptr = endptr;
	}
    }

    go = 1;
    while(go){
	begptr = line;

	tmp_nCols = 0;
	while(1){
	    val = strtod(begptr, &endptr);

	    if(val == 0 && begptr == endptr)
		break;

	    tmp_nCols++;
	    begptr = endptr;
	}
	if( tmp_nCols != nCols){
	    fprintf(stderr,"Error: Model file contains lines with different number of colums.\n");
	    exit(1);
	}

	if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) {
	    go = 0;
	}else
	    nLines++;
    }
  
    //nY = nCols;
    /* nLines = nLines; */

    /* learned weight vector */
    wLocal = (double*)calloc(nLines*nCols,sizeof(double));
    if(wLocal == NULL){
	fprintf(stderr,"Not enough memory for matrix wLocal.\n");
	exit(1);
    }


    fseek(fid,0,SEEK_SET);
    for(i=0; i < nLines; i++){
	if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL )  {
	    fprintf(stderr,"Model file corrupted.\n");
	    exit(1);
	}

	begptr = line;
	for(j=0; j < nCols; j++) {
	    val = strtod(begptr, &endptr);

	    if(val == 0 && begptr == endptr) {
		fprintf(stderr,"Model file corrupted.\n");
		exit(1);
	    }
	    begptr = endptr;
        
        
	    wLocal[LIBOCAS_INDEX(i,j,nLines)] = val;
	}
    }


    fclose(fid);
    free(line);
    toReturn->wMatrix     =wLocal;
    toReturn->nDimensions =nLines;
    toReturn->nColumns    =nCols;

    return toReturn;
}



/** 

   This function uses the model to predict the dataToPredict which has the 8 or 12 values
   for the current cycle.
   

   Again, his looks like C code because it was taken from linclass.c
   in the libocas package.
	 
   @param[in]     _dataToPredict Pointer to an array of char with the filename
   @param[in]     _nDimInModel   int with the number of dimensions in the model (should be 8 or 12)
   @param[in]     _matrixToUse   Pointer to MODELDATA struct which contains the model for the current cycle

   @return Pointer to a RETURNEDRESULTS object created on the heap with the predicted label and values

*/
RETURNEDRESULTS * singleRowData(const double * dataToPredict,const int nDimInModel,const MODELDATA * matrixToUse){

    //uint32_t i, j;
    int  i, j;
    long max_dim = 0;
    double dfce, max_dfce;
    RETURNEDRESULTS * predictionToReturn;

    predictionToReturn         = (RETURNEDRESULTS *)calloc(1, sizeof(RETURNEDRESULTS));
    predictionToReturn->scores = (double *)calloc(matrixToUse->nColumns, sizeof(double));


    if( predictionToReturn->scores == NULL ){
	fprintf(stderr,"Not enough memmory to allocate scores.\n");
	exit(1);
    }


    max_dim = LIBOCAS_MAX(max_dim,nDimInModel);
    if(matrixToUse->nDimensions != nDimInModel ){
	cerr << "The number of dimensions in the model "<<matrixToUse->nDimensions<<" does not correspond to the number specified in the argument "<<nDimInModel<<endl;
	exit(1);
    }

    max_dfce = -LIBOCAS_PLUS_INF;
    for(j=0; j < matrixToUse->nColumns; j++) {
	dfce = 0;
	for(i=0; i < nDimInModel; i++) {
	    if(i <  matrixToUse->nDimensions)
		dfce += dataToPredict[i]*matrixToUse->wMatrix[LIBOCAS_INDEX(i,j,matrixToUse->nDimensions)];
	}
	predictionToReturn->scores[j] = dfce ;

	if(max_dfce < dfce){
	    max_dfce = dfce;
	    predictionToReturn->pred_label = j+1;
	}
    }

    return predictionToReturn;
}



void read_model_index( vector<Tmodel> &models,            //to store the SVM models for the base calling
		       const char* modelindex )
{
    if(verbosity>=2) cout << "Reading model index file: " << modelindex << endl; 
    string stringPrefix1="model";
    string stringPrefix2="recal";
    string stringSuffix1=".txt";
    string stringSuffix2=".svm";
    string stringSuffix3=".par";

    ifstream ifs(modelindex);
    string line ;
    while( getline( ifs, line ) ){
	//cout<<"line "<<line <<endl;
	istringstream iss( line ) ;
	string fname, mtype ;
	int cycle ;

	if( iss >> cycle >> fname >> mtype ) {
	    cyclesInModelFile=cycle;
	    models.push_back( Tmodel() ) ;
	    Tmodel *newModel = &models.back() ;

	    if( 0 == access( fname.c_str(), R_OK ) ){
		
		newModel->matrix = readW(fname.c_str());


		if(recalibrateQuals){    
		    if(!quick)
			newModel->qualScoreSVM = new MODELDATA * [4];
		    newModel->linearRegParams  = new LinearRegParams * [4];
		    string filePrefix          = string(fname);

		    filePrefix.replace(    filePrefix.rfind(stringPrefix1),     stringPrefix2.length(),stringPrefix2); 
		    filePrefix.replace(    filePrefix.rfind(stringSuffix1),     stringSuffix2.length(),""); 
		    for(int myBaseToRead=1;myBaseToRead<=4;myBaseToRead++){
			if(!quick){
			    ostringstream fileWithSVM;
			    fileWithSVM<<filePrefix<<"_"<<myBaseToRead<<stringSuffix2;
			    // MODELDATA * tempmod=;
			    // newModel->qualScoreSVM[myBaseToRead-1]    =*( tempmod );  // we should 
			    newModel->qualScoreSVM[myBaseToRead-1]    =readW(fileWithSVM.str().c_str());
			    // free(tempmod->wMatrix);
			    // free(tempmod);


			}

			ostringstream fileWithLRP;
			fileWithLRP<<filePrefix<<"_"<<myBaseToRead<<stringSuffix3;
			newModel->linearRegParams[myBaseToRead-1] =  readLinearParams(fileWithLRP.str().c_str());	       
			
		    }
		}

	    }else {
		char *fn, *dir, *base ;
		dir = strdup( modelindex ) ;
		base = strdup( fname.c_str() ) ;
		if( -1 == asprintf( &fn, "%s/%s", dirname( dir ), basename( base ) ) )
		    throw "out of memory" ;
		newModel->matrix = readW(fn);


		if(recalibrateQuals){
		    if(!quick)
			newModel->qualScoreSVM   = new MODELDATA * [4];
		    newModel->linearRegParams= new LinearRegParams * [4];
		    string filePrefix        = string(fname);

		    filePrefix.replace(    filePrefix.rfind(stringPrefix1),     stringPrefix2.length(),stringPrefix2); 
		    filePrefix.replace(    filePrefix.rfind(stringSuffix1),     stringSuffix2.length(),""); 
		    for(int myBaseToRead=1;myBaseToRead<=4;myBaseToRead++){
			if(!quick){
			    ostringstream fileWithSVM;
			    fileWithSVM<<filePrefix<<"_"<<myBaseToRead<<stringSuffix2;
			    newModel->qualScoreSVM[myBaseToRead-1]   =readW(fileWithSVM.str().c_str());		      
			}
			ostringstream fileWithLRP;
			fileWithLRP<<filePrefix<<"_"<<myBaseToRead<<stringSuffix3;
			newModel->linearRegParams[myBaseToRead-1]=  readLinearParams(fileWithLRP.str().c_str());	       

		    }
		}



		free( fn ) ;
		free( base ) ;
		free( dir ) ;
	    }

	    newModel->type = mtype[0];

	    for( int i = 0 ; i != 4 ; ++i )
		for( int j = 0 ; j != 4 ; ++j )
		    iss >> newModel->factors[i][j] 
			>> newModel->sfactors[i][j]
			>> newModel->rclassifier[i][j] ;

	    // did we successfully parse enough factors?
	    if( iss ) 
		newModel->use_factor = true;
	    else 
		cout << "Did not have parameters for quality scores. Falling back to defaults." << endl;
	}else{
	    cerr<<"Line "<<line<<" did not parse"<<endl;
	    exit(1);
	} //end for 
    }//end for each line (cycle)

    if( verbosity >= 2 ) { cout << "Read "<<models.size()<<" Models." << endl; }
    if( models.size() == 0 ) throw "No models read" ;
    if( models.front().type != 'B' ) throw "first model must be of type 'M'" ;
    if( models.back().type != 'E' ) throw "last model must be of type 'E'" ;
}

void print_help()
{
    cout << "Usage: svm_struct_classify [options] firecrest_file model_index output_file (or directory for BCL)" << endl << endl;
    cout << "options: -h           -> this help" << endl;
    cout << "         -e EXPID     -> Name of experiment used for sequences (default TEST)"  << endl;
    cout << "         -i FILENAME  -> Index File needed when parsing IPAR instead of firecrest output"  << endl;
    cout << "         -t [1..3]    -> Type of coordinate transformation (1: round, 2: float_abs, 3: round_shift)"  << endl;
    cout << "         -o INT       -> Quality score offset (default 33)"  << endl;
    cout << "         -m INT       -> Maximum quality score reported (default "<<max_qscore<<")"  << endl;
    cout << "         -n INT       -> Minimum quality score reported (default "<<min_qscore<<")"  << endl;
    cout << "         -f INT       -> Output format (default 0: FastQ, 1: 4Q 2: BAM)"  << endl;
    cout << "         -v [0..3]    -> verbosity level (default 2)"  << endl;
    cout << "         -r           -> use recalibration (default false)"  << endl;
    cout << " Specifying cycles:    "<<endl;
    cout << "         -1           -> number of forward cycles"  << endl;
    cout << "         -2           -> number of reverse cycles"  << endl;
    cout << "         -3           -> number of index 1 cycles"  << endl;
    cout << "         -4           -> number of index 2 cycles"  << endl;

};

int main_(int argc, char**argv)
{
    const char *firecrest = 0 ;
    const char *modelindex = "svm_models" ;
    const char *outputfile = "svm_predictions" ;
    numberOfClustersFound=0;
    cyclesInModelFile=0;
    enum { i_firecrest, i_ipar, i_cif } input_mode = i_firecrest ;
    const char *IPAR_index = 0 ;
    const char* experiment = "TEST" ;
    bool isPairedEnd;
    bool hasIndex1;
    bool hasIndex2;

    recalibrateQuals=false;    
    quick           =false;    
    verbosity = 2 ;

    int forwardl=0;
    int reversel=0;
    int index1l =0;
    int index2l =0; 
    int numberOfCyclesCmdArg=0;
    int totalCycles;

    // XXX here be dragons.  How about getopt_long or popt instead?
    int i ;
    for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
	switch ((argv[i])[1]){

	case 'h': print_help(); return 0 ;
	case 'v': i++; verbosity=atol(argv[i]); break;
	case 'e': i++; experiment = argv[i] ; break;
	case 'i': i++; input_mode = i_ipar; IPAR_index = argv[i] ; break;
	case 't': i++; coordtype = atoi(argv[i]); break;
	case 'f': i++; outformat = atoi(argv[i]); break;
	case 'o': i++; qscore_offset = atoi(argv[i]); break;
	case 'm': i++; max_qscore = double(atoi(argv[i])); break;
	case 'n': i++; min_qscore = double(atoi(argv[i])); break;
	case 'r': recalibrateQuals=true; break;
	case 'q': quick=true; break;
	case '1': i++; forwardl = atoi(argv[i]); numberOfCyclesCmdArg++; break;
	case '2': i++; reversel = atoi(argv[i]); numberOfCyclesCmdArg++; break;
	case '3': i++; index1l  = atoi(argv[i]); numberOfCyclesCmdArg++; break;
	case '4': i++; index2l  = atoi(argv[i]); numberOfCyclesCmdArg++; break;

	    //case '-': parse_struct_parameters_classify(argv[i],argv[i+1]);i++; break;
	default: cout << endl << "Unrecognized option: " << argv[i] << endl;
	    print_help();
	    return 1;
	}
    }


    if((i+1)>=argc) {
	cout<<"Not enough input parameters"<<endl;
	print_help();
	return 1;
    }

    if(numberOfCyclesCmdArg != 4){
	cout<<"Must enter all cycles for every part of the sequence using -1,-2,-3,-4"<<endl<<"Use 0 if not in use e.g. -2 0 for a single end or -4 for a single index"<<numberOfCyclesCmdArg<<endl;
	print_help();
	return 1;
    }

    isPairedEnd = (reversel != 0);
    hasIndex1   = (index1l  != 0);
    hasIndex2   = (index2l  != 0);
    totalCycles=forwardl+reversel+index1l+index2l;
    if ((coordtype > 3) || (coordtype < 1)){
	cout << endl << "Coordinate type out of range: " << coordtype << endl;
	print_help();
	return 1;
    }

    firecrest = argv[i];
    modelindex = argv[i+1];
    if((i+2)<argc) {
	outputfile = argv[i+2];
    }
  
    if ((outformat < 0) || (outformat > 2)){
	cout << endl << "Output format out of range: " << outformat << endl;
	print_help();
	return 1;
    }

    vector<Tmodel>            models;
    read_model_index( models,modelindex ) ;

    double myValues[12] ;

    if( input_mode == i_ipar ) {
	struct stat the_stat ;
	if( !stat( firecrest, &the_stat ) && S_ISDIR( the_stat.st_mode ) )
	    input_mode = i_cif ;
    }

    auto_ptr<Input> input(
			input_mode == i_firecrest ? new_firecrest_input( models.size(), firecrest ) :
			input_mode == i_ipar ? new_ipar_input( models.size(), IPAR_index, coordtype, firecrest ) :
			input_mode == i_cif ? new_cif_input( models.size(), new_posn_input( IPAR_index, coordtype ), firecrest ) :
			throw "messed up input mode"
			  ) ;

    ofstream outfile;

    if(outformat != 2 ){ //not bam
	if(verbosity>=1) cout << "Opening output file: " << outputfile << endl; 
	outfile.open(outputfile,ios_base::trunc);
	if( !outfile ) throw "unable to open output file "+string(outputfile);
    }else{
	string outputfiles = string( outputfile );
	bw  = new BAMWriter ( outputfiles );
    }

    // iteration over clusters (originally one cluster == one line)
    // unsigned int numberClusters=0;

    int lane, tile, posx, posy ;
    while( input->next_cluster( &lane, &tile, &posx, &posy ) )	{
	 // numberClusters++;
	vector<Tmodel>::iterator modelitr = models.begin();
	string seq,qualB,qualA,qualC,qualG,qualT;
	double iSA, iSC, iSG, iST, iSAA, iSAC, iSAG, iSAT, iSBA=0, iSBC=0, iSBG=0, iSBT=0;

	// pre-load first cycle: the first model already needs data for
	// the second cycle
	if( !input->next_cycle( &iSA, &iSC, &iSG, &iST ) ) throw "need at least one cycle of data" ;

	// intensities for *next* cycle (from the POV of the current model)	
	while( input->next_cycle( &iSAA, &iSAC, &iSAG, &iSAT ) ) {
	    if( modelitr == models.end() ) throw "need at least as many models as cycles" ;
	    if( !iSA && !iSC && !iSG && !iST ) {
		// no light at all.  no need to pester the SVM with it.
		append_zero_quality( seq, qualB, qualA, qualC, qualG, qualT ) ;
	    } else {
		RETURNEDRESULTS * y;
		switch( modelitr->type ) 
		    {
		    case 'M':
			myValues[0]=iSA;
			myValues[1]=iSC;
			myValues[2]=iSG;
			myValues[3]=iST;
			myValues[4]=iSBA;
			myValues[5]=iSBC;
			myValues[6]=iSBG;
			myValues[7]=iSBT;
			myValues[8]=iSAA;
			myValues[9]=iSAC;
			myValues[10]=iSAG;
			myValues[11]=iSAT;
			y=singleRowData(myValues,12,modelitr->matrix);
			break;

		    case 'E':
			myValues[0]=iSA;
			myValues[1]=iSC;
			myValues[2]=iSG;
			myValues[3]=iST;
			myValues[4]=iSBA;
			myValues[5]=iSBC;
			myValues[6]=iSBG;
			myValues[7]=iSBT;
			y=singleRowData(myValues,8,modelitr->matrix);
			break;

		    case 'B':
			myValues[0]=iSA;
			myValues[1]=iSC;
			myValues[2]=iSG;
			myValues[3]=iST;
			myValues[4]=iSAA;
			myValues[5]=iSAC;
			myValues[6]=iSAG;
			myValues[7]=iSAT;
			y=singleRowData(myValues,8,modelitr->matrix);
			break;

		    default:
			throw string( "found unknown model type '" ) + modelitr->type + '\'';
		    }
		append_seq_quality(seq,qualB,qualA,qualC,qualG,qualT,y,&*modelitr/*,outfileBCL[currentCycle++]*/);

		free(y->scores);
		free(y);
	    }
	    // shift data back, we'll soon enter the next cycle
	    iSBA = iSA; iSBC = iSC; iSBG = iSG; iSBT = iST;
	    iSA = iSAA; iSC = iSAC; iSG = iSAG; iST = iSAT;
	    ++modelitr;
	}

	// special casing of last cycle
	if( modelitr == models.end() ) throw "need at least as many models as cycles" ;
	if( modelitr->type != 'E' ) throw "last model must be of type 'E'" ;

	if( !iSA && !iSC && !iSG && !iST ){
	    // no light at all.  no need to pester the SVM with it.
	    append_zero_quality( seq, qualB, qualA, qualC, qualG, qualT /*,outfileBCL[currentCycle++]*/ ) ;
	} else {
	    RETURNEDRESULTS * y;
	    myValues[0]=iSA;
	    myValues[1]=iSC;
	    myValues[2]=iSG;
	    myValues[3]=iST;
	    myValues[4]=iSBA;
	    myValues[5]=iSBC;
	    myValues[6]=iSBG;
	    myValues[7]=iSBT;
	    y=singleRowData(myValues,8,modelitr->matrix);

	    append_seq_quality(seq,qualB,qualA,qualC,qualG,qualT,y,&*modelitr /*,outfileBCL[currentCycle++]*/ ) ;

	    free(y->scores);
	    free(y);
	}
	stringstream nameseq;
	nameseq<<experiment << ':' << lane << ':' << tile << ':' << posx << ':' << posy ;
	string nameseqs=nameseq.str();

        //outfile << '@' << experiment << ':' << lane << ':' << tile << ':' << posx << ':' << posy << '\n' ;

	if(int(seq.size()) != totalCycles){
	    cerr<<"Error: sequence "<<seq<<" does not have the same number of cycles are the total specified "<<totalCycles<<endl;
	    return 1;
	}

        if( outformat == 0 ){//fastq
	    outfile << '@' <<nameseqs<<"\n";
            outfile << seq << "\n+\n" << qualB << '\n' ; 
	}else{
	    if( outformat == 1 ){ //4q
		outfile << '@' <<nameseqs<<"\n";
		outfile << "SQ " << seq << "\nQA " << qualB
			<< "\n*A " << qualA << "\n*C " << qualC << "\n*G " << qualG << "\n*T " << qualT << '\n' ;
	    }else{
		if(isPairedEnd){
		    // cerr<<"seq "<<seq<<endl;
		    // cerr<<"f   "<<seq.substr(0,forwardl) <<endl;
		    // cerr<<"r   "<<seq.substr(forwardl+index1l,reversel)<<endl;
		    // cerr<<"i1  "<<((index1l==0)?(""):(seq.substr(forwardl,index1l)))<<endl;
		    // cerr<<"i2  "<<((index2l==0)?(""):(seq.substr(forwardl+index1l+reversel,index2l)))<<endl;

		    bw->writePairedSequence(nameseqs,// name,
					    seq.substr(0,forwardl), // seq,
					    qualB.substr(0,forwardl), //qual
					    seq.substr(forwardl+index1l,reversel), // seq,
					    qualB.substr(forwardl+index1l,reversel), //qual
					    ((index1l==0)?(""):(seq.substr(forwardl,index1l))),
					    ((index1l==0)?(""):(qualB.substr(forwardl,index1l))),
					    ((index2l==0)?(""):(seq.substr(forwardl+index1l+reversel,index2l))),
					    ((index2l==0)?(""):(qualB.substr(forwardl+index1l+reversel,index2l))),
					    hasIndex1,
					    hasIndex2);
		   
		}else{
		    bw->writeSingleSequence(nameseqs,// name,
					    seq.substr(0,forwardl), // seq,
					    qualB.substr(0,forwardl), //qual
					    ((index1l==0)?(""):(seq.substr(forwardl,index1l))),
					    ((index1l==0)?(""):(qualB.substr(forwardl,index1l))),
					    ((index2l==0)?(""):(seq.substr(forwardl+index1l+reversel,index2l))),
					    ((index2l==0)?(""):(qualB.substr(forwardl+index1l+reversel,index2l))),
					    hasIndex1,
					    hasIndex2);
		}
	    }
	}
	if(outformat != 2 ) //not bam, the not open is handled in the bam writer module
	    if( !outfile ) throw "error writing output" ; 
	// if(numberClusters>100)
	//     break;
    }

    outfile.close();

    for_each( models.begin(), models.end(), FreeModel() ) ;

    if(outformat == 2 ){ //not bam
	bw->close();
	delete bw;
    }

    return 0 ;
}

int main(int argc, char**argv)
{
    try { return main_( argc, argv ) ; }
    catch( const string& e ) { cerr << e << endl ; }
    catch( const char* e ) { cerr << e << endl ; }
    catch( const exception& e ) { cerr << e.what() << endl ; }
    catch( ... ) { cerr << "Ayayay!!!" << endl ; }
    return 1 ;
}
