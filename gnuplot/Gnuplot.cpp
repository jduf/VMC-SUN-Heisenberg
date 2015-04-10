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

void Gnuplot::title(std::string const& title){plot_+="set title '"+title+"'\n";}

void Gnuplot::range(std::string const& axis, std::string const& a, std::string const& b){plot_+="set "+axis+"range ["+a+":"+b+"]\n";}
void Gnuplot::range(std::string const& axis, double const& a, double const& b){range(axis,my::tostring(a),my::tostring(b));}
void Gnuplot::range(std::string const& axis, double const& a, std::string const& b){range(axis,my::tostring(a),b);}
void Gnuplot::range(std::string const& axis, std::string const& a, double const& b){range(axis,a,my::tostring(b));}
void Gnuplot::range(std::string const& axis){plot_ += "unset "+axis+"range\n"; }

void Gnuplot::margin(std::string const& l, std::string const& r, std::string const& t, std::string const& b){
	plot_ += "set lmargin at screen "+l+"\n";
	plot_ += "set rmargin at screen "+r+"\n";
	plot_ += "set tmargin at screen "+t+"\n";
	plot_ += "set bmargin at screen "+b+"\n";
}

void Gnuplot::label(std::string const& axis, std::string const& l, std::string const& options){ plot_ += "set "+axis+"label '"+l+"' "+options+"\n"; }
void Gnuplot::label(std::string const& axis){ plot_ += "unset "+axis+"label\n"; }

void Gnuplot::operator=(std::string const& s){ plot_ = s + "\n"; }
void Gnuplot::operator+=(std::string const& s){ plot_ += s + "\n"; }

void Gnuplot::save_file(){
	IOFiles w_gp(path_ + filename_+".gp",true);
	w_gp<<plot_<<IOFiles::endl;
}

void Gnuplot::create_image(bool silent){
	std::string texfile(filename_);
	size_t pos(texfile.find("."));
	while (pos != std::string::npos) {
		texfile.replace(pos,1,"");
		pos = texfile.find(".");
	}

	Linux command;
	command(Linux::gp2latex("/tmp/"+texfile,path_,filename_));
	if(!command.status()){
		command(Linux::pdflatex("/tmp/",texfile),silent);
		command(Linux::pdf2png("/tmp/" + texfile, path_ + filename_));
		command("mv /tmp/" + texfile + ".pdf " + path_ + filename_ + ".pdf");
		command("rm /tmp/" + texfile + "*");
	} else {
		std::cerr<<"Gnuplot::create_image : can't create a plot"<<std::endl;
	}
}
