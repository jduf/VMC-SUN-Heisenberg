#ifndef DEF_FUNCTION
#define DEF_FUNCTION

#include<cmath>
#include<string> /*allow the use of sleep ?!*/
#include"Vector.hpp"

class Function{
	public:
		double f(Vector<double> x);
		static double fsol[9];
};

double Function::fsol[9]={1,-1,1.0/2.0,-1.0/6.0,1.0/24.0,-1.0/120.0,1.0/720.0,-1.0/5040,1.0/40320.0};
double Function::f(Vector<double> x){
	double d(0.0);
	for(unsigned int i(0);i<x.size();i++){
		d += (x(i)-fsol[i])*(x(i)-fsol[i]);
	}
	return d;
}
#endif
