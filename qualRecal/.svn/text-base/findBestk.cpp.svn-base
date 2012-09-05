#include <iostream>
#include <iomanip>
#include <fstream> 
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

const double dMax = std::numeric_limits<double>::max();

#include "LinRegFixed.h"


double findBestKy(vector<double> xPrek,
		  vector<double> yPrek,
		  vector<double> xPostk,
		  vector<double> yPostk,
		  double k,
		  double step,
		  double stopCriteria,
		  double lowerBoundYval,
		  double upperBoundYval,
		  bool forcePositive,
		  bool forceNegative){
    double lowerBoundY=lowerBoundYval;
    double upperBoundY=upperBoundYval;
    double overallBestkY=10.0;

    while(1){
	double incrY=(upperBoundY-lowerBoundY)/double(step);
	double kY=lowerBoundY+incrY;
	double minKy=-1;
	double minErrorY=dMax;

	for(int j=0;j<(step-2);j++){
	    LinRegFixed lrprek  (xPrek,  yPrek ,k,kY,forcePositive,forceNegative);
	    LinRegFixed lrpostk (xPostk, yPostk,k,kY,forcePositive,forceNegative);
	    double totalSQ =lrprek.getSquareSum() + lrpostk.getSquareSum();
	    if(totalSQ<minErrorY){
		minErrorY=totalSQ;
		minKy=kY;
	    }
	    kY+=incrY;
	}//end for each step j


	lowerBoundY=minKy-incrY;
	upperBoundY=minKy+incrY;
	if(abs(upperBoundY-lowerBoundY) < stopCriteria){
	    overallBestkY=minKy;
	    break;
	}

    }
    return overallBestkY;
}

int main(int argc,char * argv[]){
    vector<double> x;
    vector<double> y;
    bool forcePositive=false;
    bool forceNegative=false;

    if(argc==0){
	cout<<"Usage:"<<endl
	    <<"./findBestk [file] [flag]"<<endl
	    <<"where:"<<endl
	    <<"[flag] = 0 to force positive slopes"<<endl
	    <<"       = 1 to force negative slopes"<<endl;
	exit(1);
    }
    string fileName1=argv[1];
    int flag        =atoi(argv[2]);
    if(flag == 0){
	forcePositive=true;
    }else{
	if(flag == 1){
	    forceNegative=true;
	}else{
	    cout<<"wrong flag value"<<endl
		<<"Usage:"<<endl
		<<"./findBestk [file] [flag]"<<endl
		<<"where:"<<endl
		<<"[flag] = 0 to force positive slopes"<<endl
		<<"       = 1 to force negative slopes"<<endl;
	    exit(1);
	}
    }

    int lineCounter=0;
    ifstream valuesin;
    string line;
    // double yMax;

    valuesin.open(fileName1.c_str(), ios::in);    // open the streams
    if (valuesin.is_open()) {
	while ( getline (valuesin,line) ){
	    lineCounter++;
	    istringstream stringStr( line ) ;
	    double xval;
	    double yval;
	    double errorR;
	    double yvalRound;

	    if( stringStr >> xval >> errorR >> yvalRound >> yval ) {
		// if(x.size() == 0)
		//     yMax=yval;
		// else
		//     if(yval>yMax)
		// 	yMax=yval;

		x.push_back(xval);
		y.push_back(yval);
	    }else{
		cerr<<"Verify the format of the file "<<fileName1<<" wrong line #"<<lineCounter<<" line=#"<<line<<"#"<<endl;
		return 1;
	    }   
	}
	valuesin.close();
    }else 
	cout << "Unable to open file "<<fileName1<<endl;







    double lowerBound=*min_element(x.begin(),x.end());
    double upperBound=*max_element(x.begin(),x.end());
    double lowerBoundYval=*min_element(y.begin(),y.end());
    double upperBoundYval=*max_element(y.begin(),y.end());

    int step=10;
    double stopCriteria=0.1;
    double overallBestk=10.0;

    while(1){
	double incr=(upperBound-lowerBound)/double(step);
	double k=lowerBound+incr;
	double minError=dMax;
	double minK=-1;

	for(int i=0;i<(step-2);i++){
	    vector<double> xPrek;
	    vector<double> yPrek;
	    vector<double> xPostk;
	    vector<double> yPostk;

	    for(int j=0;j<x.size();j++){
	    	if(x[j] <=k){
	    	    xPrek.push_back(x[j]);
	    	    yPrek.push_back(y[j]);
	    	}else{
		    xPostk.push_back(x[j]);
		    yPostk.push_back(y[j]);

		}
	    }


	    if(xPrek.size()  != 0 &&
	       xPostk.size() != 0 ){

		double overallBestkY=findBestKy(xPrek,yPrek,xPostk,yPostk,k,step,stopCriteria,lowerBoundYval,upperBoundYval,forcePositive,forceNegative);
		LinRegFixed lrprek  (xPrek,  yPrek  ,k,overallBestkY,forcePositive,forceNegative);
		LinRegFixed lrpostk (xPostk, yPostk ,k,overallBestkY,forcePositive,forceNegative);

		double totalSQ =lrprek.getSquareSum() + lrpostk.getSquareSum();

		if(totalSQ<minError){
		    minError=totalSQ;
		    minK=k;
		}
	    }
	    
	    k+=incr;
	} //end for each step i

	lowerBound=minK-incr;
	upperBound=minK+incr;
	if(abs(upperBound-lowerBound) < stopCriteria){
	    overallBestk=minK;
	    break;
	}

    }

    vector<double> xPrek;
    vector<double> yPrek;
    vector<double> xPostk;
    vector<double> yPostk;

    for(int j=0;j<x.size();j++){
	if(x[j] <=overallBestk){
	    xPrek.push_back(x[j]);
	    yPrek.push_back(y[j]);
	}else{
	    xPostk.push_back(x[j]);
	    yPostk.push_back(y[j]);   
	}
    }


    double overallBestkY=findBestKy(xPrek,yPrek,xPostk,yPostk,overallBestk,step,stopCriteria,lowerBoundYval,upperBoundYval,forcePositive,forceNegative);

    LinRegFixed lrprek  (xPrek,  yPrek, overallBestk,overallBestkY,forcePositive,forceNegative);
    LinRegFixed lrpostk (xPostk, yPostk,overallBestk,overallBestkY,forcePositive,forceNegative);

    cout<<overallBestk<<endl;
    cout << lrprek.getA() << endl;
    cout << lrprek.getB() << endl;
    cout << lrpostk.getA() << endl;
    cout << lrpostk.getB() << endl;


    return 0;
}
