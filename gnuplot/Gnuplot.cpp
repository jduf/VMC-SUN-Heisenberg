#include "Gnuplot.hpp"

Gnuplot::Gnuplot(std::string filename, std::string type):
	filename_(filename),
	plottype_(type),
	preplot_(""),
	plot_("")
{ }

Gnuplot::~Gnuplot(){
	Write w_gp(filename_+".gp");
	w_gp<<preplot_<<Write::endl
		<<plottype_<<" "<<plot_<<Write::endl;
}

void Gnuplot::save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Vector<double> const& z){
	Matrix<double> data(x.size(),3);
	for(unsigned int i(0);i<x.size();i++){
		data(i,0) = x(i);
		data(i,1) = y(i);
		data(i,2) = z(i);
	}

	plot_ += "'" + data_file + "'";
	Write w_data(data_file);
	w_data<<data;
}

void Gnuplot::save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Matrix<double> const& z){
	plot_ += "'" + data_file + "'";
	Write w_data(data_file);
	for(unsigned int i(0);i<x.size();i++){
		for(unsigned int j(0);j<y.size();j++){
			w_data<<x(i)<<" "<<y(j)<<" "<<z(i,j)<<Write::endl;
		}
		w_data<<Write::endl;
	}
}

void Gnuplot::save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y){
	Matrix<double> data(x.size(),2);
	for(unsigned int i(0);i<x.size();i++){
		data(i,0) = x(i);
		data(i,1) = y(i);
	}

	plot_ += "'" + data_file + "'";
	Write w_data(data_file);
	w_data<<data;
}

void Gnuplot::save_data(std::string data_file, Matrix<double> const& z){
	plot_ += "'" + data_file + "'";
	Write w_data(data_file);
	w_data<<z;
}

void Gnuplot::add_plot_param(std::string s){
	plot_ += " " + s;
}

void Gnuplot::preplot(std::string s){ 
	preplot_ += s + "\n";
}


