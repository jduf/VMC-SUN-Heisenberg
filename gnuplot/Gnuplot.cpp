#include "Gnuplot.hpp"

Gnuplot::Gnuplot(std::string const& path, std::string const& filename):
	path_(path),
	filename_(filename),
	plot_("")
{
	if(path_[0] != '/'){
		Linux c;
		path_ = c.pwd() + path_;
	}
}

Gnuplot::~Gnuplot(){}

void Gnuplot::save_file(){
	IOFiles w_gp(path_ + filename_+".gp",true);
	w_gp<<plot_<<IOFiles::endl;
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
