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

void PSTricks::put(double x, double y, std::string s, std::string options){
	s_ += "\\rput" + options + "("+tostring(x)+","+tostring(y)+"){"+s+"}\n";
}

void PSTricks::pie(Vector<double> const& x, double r, std::string options){
	s_ += "\\psChart["+options+",linestyle=none]{";
	for(unsigned int i(0); i<x.size(); i++){
		s_ += tostring(x(i));
		if(i+1<x.size()){ s_ += ","; }
	}
	s_ += "}{}{"+tostring(r)+"}\n";
}

void PSTricks::save(bool silent, bool pdf){
	{/*to make sure that the file w is closed after the brackets*/
		IOFiles w(path_+filename_ + ".tex",true);
		s_ += "\\end{document}\n";
		w<<s_;
	}

	Linux command;
	if(silent){ command("latex --shell-escape " + filename_ + ".tex > /dev/null 2> /dev/null");}
	else{ command("latex --shell-escape " + filename_ + ".tex");}
	if(pdf){ command("dvipdf " + filename_ + ".dvi " + path_ + filename_ + ".pdf "); }
	command("mv " + filename_ + "-1.png " + path_ + filename_ + ".png");
	command("rm *.dvi *.aux *.log *.ps" );
}

