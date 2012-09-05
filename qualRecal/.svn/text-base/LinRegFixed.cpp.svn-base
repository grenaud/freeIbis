#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "LinRegFixed.h"

LinRegFixed::LinRegFixed(vector<double> x, vector<double> y,double xk,double yk,bool forcePositive,bool forceNegative,bool eliminateBegin){

    if(x.size() != y.size() ){
	cout<<"The size of the x vector ("<<x.size() <<") is not the same as the y vector ("<<y.size() <<")"<<endl;
	exit(1);
    }

    if(forcePositive == forceNegative && 
       forcePositive){
	cout<<"Cannot force the slope to be both  negative and positive"<<endl;
	exit(1);
    }

    double size=double(x.size());
    double numerator=0.0;
    double denominator=0.0;
    for(int i=0;i<x.size();i++){
	numerator+= ( (yk-y[i])*(xk-x[i]) );
	denominator += pow(xk-x[i],2.0);
    }

    b=numerator/denominator;
    a=yk-b*xk;

    if(forcePositive){
	if(b<0){
	    a=yk;
	    b=0.0;
	}
    }

    if(forceNegative){
	if(b>0){
	    a=yk;
	    b=0.0;
	}
    }

    vector<double> xNoBegin;
    vector<double> yNoBegin;
    double fractionOfBeginning  = 0.25;
    double maxvalForLeastSquare = 20.0;

    sumSquares=0;
    for(int i=0;i<x.size();i++){
	double lqtoadd=pow( (y[i]-(a+b*x[i])) ,2.0);
	sumSquares+=lqtoadd;
	if(eliminateBegin){
	    if(i < int(0.25*x.size())){
		if(lqtoadd < maxvalForLeastSquare){
		    xNoBegin.push_back(x[i]);
		    yNoBegin.push_back(y[i]);
		}
	    }else{
		xNoBegin.push_back(x[i]);
		yNoBegin.push_back(y[i]);
	    }
	}
    }

    if(eliminateBegin){
	size=double(xNoBegin.size());
	numerator=0.0;
	denominator=0.0;
	for(int i=0;i<xNoBegin.size();i++){
	    numerator+= ( (yk-yNoBegin[i])*(xk-xNoBegin[i]) );
	    denominator += pow(xk-xNoBegin[i],2.0);
	}

	b=numerator/denominator;
	a=yk-b*xk;

	if(forcePositive){
	    if(b<0){
		a=yk;
		b=0.0;
	    }
	}

	if(forceNegative){
	    if(b>0){
		a=yk;
		b=0.0;
	    }
	}
    }

}

// LinRegFixed::~LinRegFixed(){

// }


