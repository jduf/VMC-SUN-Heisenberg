#include "Gnuplot.hpp"

Gnuplot::Gnuplot(std::string path, std::string filename, std::string type):
	path_(path),
	filename_(filename),
	plottype_(type),
	preplot_(""),
	plot_("")
{}

Gnuplot::Gnuplot(std::string path, std::string filename):
	path_(path),
	filename_(filename),
	plottype_(""),
	preplot_(""),
	plot_("")
{}

Gnuplot::~Gnuplot(){}

void Gnuplot::save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Vector<double> const& z){
	Matrix<double> data(x.size(),3);
	for(unsigned int i(0);i<x.size();i++){
		data(i,0) = x(i);
		data(i,1) = y(i);
		data(i,2) = z(i);
	}

	plot_ += "'" + data_file + "'";
	Write w_data(path_+data_file);
	w_data<<data;
}

void Gnuplot::save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Matrix<double> const& z){
	plot_ += "'" + data_file + "'";
	Write w_data(path_+data_file);
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
	Write w_data(path_+data_file);
	w_data<<data;
}

void Gnuplot::save_data(std::string data_file, Matrix<double> const& z){
	plot_ += "'" + data_file + "'";
	Write w_data(path_+data_file);
	w_data<<z;
}

void Gnuplot::add(std::string s){
	plot_ += " " + s;
}

void Gnuplot::preplot(std::string s){ 
	preplot_ += s + "\n";
}

void Gnuplot::save_file(){
	Write w_gp(path_ + filename_+".gp");
	w_gp<<preplot_<<Write::endl
		<<plottype_<<" "<<plot_<<Write::endl;
}

void Gnuplot::create_image(bool silent){
	std::string tmp(filename_);
	size_t pos(tmp.find("."));
	while (pos != std::string::npos) {
		tmp.replace(pos,1,"");
		pos = tmp.find(".");
	}

	Linux command;
	command("cd " + path_ + "; gnuplot -e \"set terminal epslatex color size 12.5cm,7.73 standalone lw 2 header \'\\\\usepackage{amsmath}\'; set output \'/tmp/" + tmp + ".tex\'\" " + path_ + filename_ + ".gp");
	if(!command.status()){
		if(silent){ command("cd /tmp/; pdflatex -shell-escape " + tmp + ".tex > /dev/null 2> /dev/null");}
		else{ command("cd /tmp/; pdflatex -shell-escape " + tmp + ".tex");}
		command("cd /tmp/; convert -density 500 -resize 20% " + tmp + ".pdf " + path_ + filename_ + ".png");
		command("mv /tmp/" + tmp + ".pdf " + path_ + filename_ + ".pdf");
		command("rm /tmp/" + tmp + "*");
	} else {
		std::cerr<<"Gnuplot::create_image : can't create a plot"<<std::endl;
	}
}

