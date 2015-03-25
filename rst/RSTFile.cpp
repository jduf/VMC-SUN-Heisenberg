#include "RSTFile.hpp"

RSTFile::RSTFile(std::string const& path, std::string const& filename):
	path_(path),
	filename_(filename)
{}

void RSTFile::save(bool pdf){
	IOFiles w(path_ + filename_ + ".rst",true);
	rst_ += RST_np_;
	w<<rst_;
	Linux command;
	command(Linux::rst2html(path_,filename_));
	if(command.status()){
		std::cerr<<"RSTFile::~RSTFile() : the command function called returns an error for the html creation"<<std::endl; 
	} else {
		if(pdf){
			command(Linux::rst2latex(path_,filename_));
			command(Linux::pdflatex(path_,filename_));
			if(command.status()){std::cerr<<"RSTFile::~RSTFile() : the command function called returns an error for the pdf creation"<<std::endl; }
		}
	}
}
