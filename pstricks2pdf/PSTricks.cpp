#include "PSTricks.hpp"

PSTricks::PSTricks(std::string path, std::string filename):
	path_(path),
	filename_(filename),
	s_(""),
	begin_end_(false)
{
	s_ +="\\documentclass[crop,pstricks,png]{standalone}\n" ;
	s_ +="\\usepackage{pstricks-add}\n";
	s_ +="\\begin{document}\n";
}

void PSTricks::begin(double const& xbl, double const& ybl, double const& xtr, double const& ytr, std::string const& imagename){
	s_ += "\\begin{pspicture}("+my::tostring(xbl)+","+my::tostring(ybl)+")("+my::tostring(xtr)+","+my::tostring(ytr)+")%"+imagename+"\n";
	begin_end_ = true;
}

void PSTricks::add(std::string const& s){
	s_ += s + "\n";
}

void PSTricks::line(std::string const& linetype, double const& x0, double const& y0, double const& x1, double const& y1, std::string const& options){
	s_ += "\\psline[" + options + "]{" + linetype + "}(" + my::tostring(x0) + "," + my::tostring(y0) + ")(" + my::tostring(x1) + "," + my::tostring(y1) + ")\n";
}

void PSTricks::polygon(Matrix<double> const& xy, std::string const& options){
	s_ += "\\pspolygon[" + options + "]";
	for(unsigned int i(0);i<xy.row();i++){
		s_ += "(" + my::tostring(xy(i,0)) + "," + my::tostring(xy(i,1)) + ")";
	}
	s_ += "\n";
}

void PSTricks::frame(double const& x0, double const& y0, double const& x1, double const& y1, std::string const& options){
	s_ += "\\psframe[" + options + "](" + my::tostring(x0) + "," + my::tostring(y0) + ")(" + my::tostring(x1) + "," + my::tostring(y1) + ")\n";
}

void PSTricks::put(double const& x, double const& y, std::string const& s, std::string const& options){
	s_ += "\\rput" + options + "("+my::tostring(x)+","+my::tostring(y)+"){"+s+"}\n";
}

void PSTricks::pie(Vector<double> const& x, double const& r, std::string const& options){
	s_ += "\\psChart["+options+",linestyle=none]{";
	for(unsigned int i(0); i<x.size(); i++){
		s_ += my::tostring(x(i));
		if(i+1<x.size()){ s_ += ","; }
	}
	s_ += "}{}{"+my::tostring(r)+"}\n";
}

void PSTricks::circle(Vector<double> const& x, double const& r, std::string const& options){
	s_ += "\\pscircle["+options+"]("+my::tostring(x(0))+","+my::tostring(x(1))+"){"+my::tostring(r)+"}\n";
}

void PSTricks::end(bool const& silent, bool const& pdf, bool const& crop){
	if(begin_end_){ s_ += "\\end{pspicture}\n"; }
	s_ += "\\end{document}\n";
	{/*to make sure that the file w is closed after the brackets*/
		IOFiles w(path_+filename_ + ".tex",true);
		w<<s_;
	}

	Linux command;
	command(Linux::latex(path_,filename_),silent);
	if(!command.status()){
		if(pdf){
			command(Linux::dvipdf(path_,filename_),silent);
			if(crop){ command(Linux::pdfcrop(path_,filename_),silent); }
			command(Linux::pdf2png(path_+filename_,path_+filename_),silent);
			command("rm "+path_+filename_+".dvi "+path_+filename_+".aux "+path_+filename_+".log "+path_+filename_+".ps "+path_+filename_+"-1.png",true);
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : Linux::latex(path_,filename_) returned an error ("<<command.status()<<")"<<std::endl; }
}
