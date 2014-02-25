#include "PSTricks.hpp"

PSTricks::PSTricks(std::string filename, std::string path):
	path_(path),
	filename_(filename),
	s_("")
{
	s_ +="\\documentclass{article}\n" ;
	s_ +="\\usepackage[pdf]{pstricks}\n";
	s_ +="\\usepackage{pstricks-add}\n";
	s_ +="\\begin{document}\n";
	s_ +="\\titlepage\n"; /*if i can find a better way. but \tiny has to work*/
}

PSTricks::~PSTricks(){
	Write w(path_ + filename_ + ".tex");
	s_ += "\\end{document}\n";

	w<<s_;

	Linux command;
	command("latex --output-directory=/tmp " + path_ + filename_ + ".tex");
	command("dvips /tmp/" + filename_ + ".dvi -o /tmp/" + filename_ + ".ps" );
	command("ps2pdf /tmp/" + filename_ + ".ps /tmp/" + filename_ + ".pdf");
	command("pdfcrop /tmp/" + filename_ + ".pdf " + path_ + filename_ + ".pdf" );
}

void PSTricks::add(std::string s){
	s_ += s + "\n";
}

void PSTricks::line(std::string linetype, double x0, double y0, double x1, double y1, std::string options){
	s_ += "\\psline[" + options + "]{" + linetype + "}(" + tostring(x0) + "," + tostring(y0) + ")(" + tostring(x1) + "," + tostring(y1) + ")\n";
}

void PSTricks::polygon(Matrix<double> const& xy, std::string options){
	s_ += "\\pspolygon[" + options + "]";
	for(unsigned int i(0);i<xy.row();i++){
		s_ += "(" + tostring(xy(i,0)) + "," + tostring(xy(i,1)) + ")";
	}
}

void PSTricks::frame(double x0, double y0, double x1, double y1, std::string options){
	s_ += "\\psframe[" + options + "](" + tostring(x0) + "," + tostring(y0) + ")(" + tostring(x1) + "," + tostring(y1) + ")\n";
}

void PSTricks::put(double x, double y, std::string s){
	s_ += "\\rput("+tostring(x)+","+tostring(y)+"){"+s+"}\n";
}

void PSTricks::pie(Vector<double> const& x, double r, std::string options){
	s_ += "\\psChart["+options+",linestyle=none]{";
	for(unsigned int i(0); i<x.size(); i++){
		s_ += tostring(x(i));
		if(i+1<x.size()){ s_ += ","; }
	}
	s_ += "}{}{"+tostring(r)+"}\n";
}
