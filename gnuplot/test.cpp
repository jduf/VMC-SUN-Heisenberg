#include "Gnuplot.hpp"

#include<iostream>

int main(){
	Vector<double> x(100,0.0);
	Vector<double> y(300,0.0);
	Matrix<double> m(x.size(),y.size(),0.0);
	for(unsigned int i(0);i<x.size();i++){
		x(i)=i*2*M_PI/x.size();
	}
	for(unsigned int i(0);i<y.size();i++){
		y(i)=i*2*M_PI/y.size();
	}
	for(unsigned int i(0);i<x.size();i++){
		for(unsigned int j(0);j<y.size();j++){
			m(i,j) = sin(x(i))+cos(y(j));
		}
	}

	Linux command;

	Gnuplot gp(command.pwd(),"bla","plot");
	gp.save_data("data.dat",x,y,m);
	gp.save_file();
	gp.create_image(false);
}
