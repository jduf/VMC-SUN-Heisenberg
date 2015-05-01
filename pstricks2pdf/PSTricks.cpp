#include "PSTricks.hpp"

PSTricks::PSTricks(std::string path, std::string filename):
	path_(path),
	filename_(filename),
	s_("")
{
	s_ +="\\documentclass[crop,pstricks,png]{standalone}\n" ;
	s_ +="\\usepackage{pstricks-add}\n";
	s_ +="\\begin{document}\n";
}

void PSTricks::add(std::string s){
	s_ += s + "\n";
}

void PSTricks::line(std::string linetype, double x0, double y0, double x1, double y1, std::string options){
	s_ += "\\psline[" + options + "]{" + linetype + "}(" + my::tostring(x0) + "," + my::tostring(y0) + ")(" + my::tostring(x1) + "," + my::tostring(y1) + ")\n";
}

void PSTricks::polygon(Matrix<double> const& xy, std::string options){
	s_ += "\\pspolygon[" + options + "]";
	for(unsigned int i(0);i<xy.row();i++){
		s_ += "(" + my::tostring(xy(i,0)) + "," + my::tostring(xy(i,1)) + ")";
	}
	s_ += "\n";
}

void PSTricks::frame(double x0, double y0, double x1, double y1, std::string options){
	s_ += "\\psframe[" + options + "](" + my::tostring(x0) + "," + my::tostring(y0) + ")(" + my::tostring(x1) + "," + my::tostring(y1) + ")\n";
}

void PSTricks::put(double x, double y, std::string s, std::string options){
	s_ += "\\rput" + options + "("+my::tostring(x)+","+my::tostring(y)+"){"+s+"}\n";
}

void PSTricks::pie(Vector<double> const& x, double r, std::string options){
	s_ += "\\psChart["+options+",linestyle=none]{";
	for(unsigned int i(0); i<x.size(); i++){
		s_ += my::tostring(x(i));
		if(i+1<x.size()){ s_ += ","; }
	}
	s_ += "}{}{"+my::tostring(r)+"}\n";
}

void PSTricks::circle(Vector<double> const& x, double r, std::string options){
	s_ += "\\pscircle["+options+"]("+my::tostring(x(0))+","+my::tostring(x(1))+"){"+my::tostring(r)+"}\n";
}

void PSTricks::save(bool silent, bool pdf, bool crop){
	{/*to make sure that the file w is closed after the brackets*/
		IOFiles w(path_+filename_ + ".tex",true);
		s_ += "\\end{document}\n";
		w<<s_;
	}

	Linux command;
	command(Linux::latex(filename_),silent);
	if(pdf){ command(Linux::dvipdf(path_,filename_),silent);
		if(crop){ command(Linux::pdfcrop(path_,filename_),silent); }
	}
	command("mv " + filename_ + "-1.png " + path_ + filename_ + ".png",silent);
	command("rm *.dvi *.aux *.log *.ps" ,silent);
}
