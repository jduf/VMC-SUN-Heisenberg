#include "Gnuplot.hpp"

#include<iostream>

int main(){
	Matrix<double> x(50,1,0.0);
	Matrix<double> y(30,1,0.0);
	Matrix<double> m(x.total(),y.total(),0.0);
	for(unsigned int i(0);i<x.total();i++){
		x(i)=i*2*M_PI/x.total();
	}
	for(unsigned int i(0);i<y.total();i++){
		y(i)=i*2*M_PI/y.total();
	}
	for(unsigned int i(0);i<x.total();i++){
		for(unsigned int j(0);j<y.total();j++){
			m(i,j) = sin(x(i))+cos(y(j));
		}
	}


	Gnuplot gp("bla","3D");
	gp.data(x,y,m);
	gp.save();
}
