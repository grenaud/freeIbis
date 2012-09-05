#ifndef LinRegFixed_h
#define LinRegFixed_h

#include <vector>

using namespace std;

class LinRegFixed{
private:
    double a;
    double b;
    double sumSquares;

public:
    LinRegFixed(vector<double> x, vector<double> y,double xk,double yk,bool forcePositive,bool forceNegative,bool eliminateBegin=false);
    //~LinRegFixed();
    double getA() const {return a;}; //get the intercept
    double getB() const {return b;}; //get the slope
    double getSquareSum() const { return sumSquares; }; //get the sqError

};
#endif
