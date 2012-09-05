#include <iostream>
#include <fstream> 
#include <iomanip> 
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;


int main (int argc, char *argv[]) {
    ifstream scoresin;
    int lineCounter=0;
    string line;
    vector<float> diffDist;
    vector<int>   mmvec;
    int requiredMismatches=200;
    int bound1  =atoi(argv[1]);
    int bound2  =atoi(argv[2]);
    int errorMax=atoi(argv[3]);
    string fileName=argv[4];

    int windowSize;

    windowSize=round(pow(10.0,errorMax/10.0));

    scoresin.open(fileName.c_str(), ios::in);    // open the streams
    if (scoresin.is_open()) {
	while ( getline (scoresin,line) ){
	    lineCounter++;
	    istringstream stringStr( line ) ;
	    float diff;
	    int    mm;
	    if( stringStr >> diff >> mm ) {
		diffDist.push_back(diff);
		mmvec.push_back(mm);
	    }else{
		cerr<<"Verify the format of the file "<<fileName<<" wrong line #"<<lineCounter<<" line=#"<<line<<"#"<<endl;
		return 1;
	    }   
	}
	scoresin.close();
    }else{
	cout << "Unable to open file "<<fileName<<endl;
    }
   
    vector<double> estimatedDistance;
    vector<double> estimatedError;
    int numberOfIterations=0;

    while (1) {
	int i=0;
	while (i<diffDist.size()) {
	   int j=0;
	   double sum=0.0;
	   int counter=0;
	   int counterMM=0;
	   bool correctWindow=false; 	   
	   while ( (i+j) <= diffDist.size() ) {
	       if (mmvec[i+j] == 0) {
		   counterMM++;
	       }
	       sum+=diffDist[i+j];

	       if (counterMM >= requiredMismatches) {
		   correctWindow=true;
		   break;
	       }
	       
	       j++;
	       counter++; 
	   }
	   
	   double smallWindowValue=(sum/double(counter));
	   double smallWindowEst  =double(double(counterMM)/double(counter));
	   if(correctWindow){
	       estimatedDistance.push_back(smallWindowValue);
	       estimatedError.push_back(smallWindowEst);
	   }else{
	       break;
	   }

	   if(int(counter/4) > 1)
	       i+=int(counter/4);
	   else
	       i+=2;

	}

	// cout<<"requiredMismatches="<<requiredMismatches<<endl;
	// cout<<"size              ="<<estimatedDistance.size() <<endl;

	if(estimatedDistance.size() < bound1 ){
	    if(requiredMismatches < 3){
		requiredMismatches--;
	    }else{
		double factor=(double(rand()%75)/100.0+0.5);
		requiredMismatches=int(requiredMismatches*factor);
	    }
	    if(requiredMismatches <= 1){
		estimatedError.clear();
		estimatedDistance.clear();
		break;		
	    }
	    estimatedError.clear();
	    estimatedDistance.clear();

	}else{
	    if(estimatedDistance.size() > bound2 ){
		if(requiredMismatches < 3){
		    requiredMismatches++;
		}else{
		    double factor=(double(rand()%100)/100.0 +1.25);
		    requiredMismatches=int(requiredMismatches*factor);
		}
		estimatedError.clear();
		estimatedDistance.clear();

	    }else{
		estimatedError.clear();
		estimatedDistance.clear();
		break;
	    }
	}
	numberOfIterations++;
	if(numberOfIterations>5000)
	    break;
    }// while(1)
    


    int i=0;
    int numberOftimesWithMaxScoreMAX=10;
    int numberOftimesWithMaxScore=0;

    while (i<diffDist.size()) {
	int j=0;
	double sum=0.0;
	int counter=0;
	int counterMM=0;
	bool correctWindow=false; 	   
	while ( (i+j) <= diffDist.size() ) {
	    if (mmvec[i+j] == 0) {
		counterMM++;
	    }
	    sum+=diffDist[i+j];
	    j++;
	    counter++; 
	    if (counterMM >= requiredMismatches) {
		correctWindow=true;
		break;
	    }
	       
	  
	    if(j>=windowSize){
		if(counterMM <= 1){
		    numberOftimesWithMaxScore++;
		    counterMM=1;
		}
		correctWindow=true;
		break;
	    }

	}
	// cout<<endl<<"sum "<<sum<<endl;
	// cout<<"counter "<<counter<<endl;
	// cout<<"counterMM "<<counterMM<<endl;

	double smallWindowValue;
	double smallWindowEst;
	
	// if(numberOftimesWithMaxScore > numberOftimesWithMaxScoreMAX){
	//     smallWindowEst=double(double(1)/double(windowSize));
	// }else{
	smallWindowValue=(sum/double(counter));
	smallWindowEst=double(double(counterMM)/double(counter));
	// }
	// cout<<"smallWindowValue "<<smallWindowValue<<endl;
	// cout<<"smallWindowEst "<<smallWindowEst<<endl;


	if(correctWindow){
	    if( !(counterMM == 1 && counter == 1) ){
		estimatedDistance.push_back(smallWindowValue);
		estimatedError.push_back(smallWindowEst);
	    }
	}else{
	    break;
	}

	if(int(counter/4) > 1)
	    i+=int(counter/4);
	else
	    i+=2;

    }


    for ( int k=0;k<estimatedDistance.size();k++){
	int estimatedInt;
	double estimatedFl;

	if(estimatedError[k] == 0){
	    estimatedInt=round(errorMax);
	    estimatedFl=double(errorMax);
	}else{
	    estimatedInt=round(-10.0*log10(estimatedError[k]));
	    estimatedFl =(-10.0*log10(estimatedError[k]));
	}
	cout  << estimatedDistance[k]<<"\t"<<estimatedError[k]<<"\t"<<estimatedInt<<"\t"<<estimatedFl<<endl;
    }

    return 0;
}

