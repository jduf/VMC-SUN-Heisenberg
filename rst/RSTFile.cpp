#include "RSTFile.hpp"

RSTFile::RSTFile(std::string const& path, std::string const& filename):
	path_(path),
	filename_(filename)
{}

void RSTFile::save(bool const& pdf, bool const& silent){
	{
		IOFiles w(path_+filename_+".rst",true,false);
		w<<rst_;
	}
	Linux command;
	command(Linux::rst2html(path_,filename_),silent);
	if(command.status()){
		std::cerr<<__PRETTY_FUNCTION__<<" : Linux::rst2html(path_,filename_) returned an error ("<<command.status()<<")"<<std::endl;
	} else {
		if(pdf){
			command(Linux::rst2latex("/tmp/"+filename_,path_,filename_),silent);
			command(Linux::pdflatex("/tmp/",filename_),silent);
			if(command.status()){ std::cerr<<__PRETTY_FUNCTION__<<" : Linux::pdf2latex(path_,filename_) returned an error ("<<command.status()<<")"<<std::endl; }
		}
	}
}
