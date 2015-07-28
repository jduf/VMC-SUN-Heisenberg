#include "Gnuplot.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

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

	Gnuplot gp(command.pwd(),"bla");
	IOFiles w("data.dat",true);
	for(unsigned int i(0);i<x.size();i++){
		for(unsigned int j(0);j<y.size();j++){
			w<<x(i)<<" "<<y(j)<<" "<<m(i,j)<<IOFiles::endl;
		}
		w<<IOFiles::endl;
	}
	gp += "set view 80,40";
	gp += "splot 'data.dat'";
	gp.save_file();
	gp.create_image(false,true);
}
