#include "Gnuplot.hpp"

Gnuplot::Gnuplot(std::string filename, std::string type):
	filename(filename),
	gp_code("")
{
	if(type == "3D"){
		gp_code = "splot '" + filename + ".dat' ";
	}
}

Gnuplot::~Gnuplot(){}

void Gnuplot::data(Matrix<double> const& x, Matrix<double> const& y, Matrix<double> const& z){
	this->gp_data = Matrix<double> (x.total()*y.total(),3);
	for(unsigned int i(0);i<x.total();i++){
		for(unsigned int j(0);j<y.total();j++){
			this->gp_data(j+i*y.total(),0) = x(i);
			this->gp_data(j+i*y.total(),1) = y(j);
			this->gp_data(j+i*y.total(),2) = z(i,j);
		}
	}
}

void Gnuplot::code(std::string s){
	this->gp_code += s + "\n";
}

void Gnuplot::save(){
	Write w_data(this->filename+".dat");
	w_data<<this->gp_data;
	Write w_gp(this->filename+".gp");
	w_gp<<this->gp_code;
}

void Gnuplot::test(){ }


