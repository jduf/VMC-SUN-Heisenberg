#include "Gnuplot.hpp"

Gnuplot::Gnuplot(std::string filename, std::string type):
	filename(filename),
	gp_code("")
{
	if(type == "3D"){
		gp_code = "splot ";
	}
	if(type == "1D"){
		gp_code = "plot ";
	}
}

Gnuplot::~Gnuplot(){}

void Gnuplot::save_data(std::string data_file, Matrix<double> const& x, Matrix<double> const& y, Matrix<double> const& z){
	Matrix<double> data;
	if( x.row() == z.total() && y.row() == z.total()){
		data.set(x.row(),3);
		for(unsigned int i(0);i<x.row();i++){
			data(i,0) = x(i);
			data(i,1) = y(i);
			data(i,2) = z(i);
			std::cout<<i<<" "<<x(i)<<" "<<data(i,0)<<std::endl;
		}
	} else {
		data.set(x.row()*y.row(),3);
		for(unsigned int i(0);i<x.row();i++){
			for(unsigned int j(0);j<y.row();j++){
				data(j+i*y.row(),0) = x(i);
				data(j+i*y.row(),1) = y(j);
				data(j+i*y.row(),2) = z(i,j);
			}
		}
	}

	this->gp_code += "'" + data_file + "'";
	Write w_data(data_file);
	w_data<<data;
}

void Gnuplot::save_data(std::string data_file, Matrix<double> const& x, Matrix<double> const& y){
	this->gp_code += data_file;

	Matrix<double> data(x.row(),2);
	for(unsigned int i(0);i<x.total();i++){
		data(i,0) = x(i);
		data(i,1) = y(i);
	}

	this->gp_code += "'" + data_file + "'";
	Write w_data(data_file);
	w_data<<data;
}

void Gnuplot::code(std::string s){
	this->gp_code += s;
}

void Gnuplot::save_code(){
	Write w_gp(this->filename+".gp");
	w_gp<<this->gp_code;
}

void Gnuplot::test(){ }


