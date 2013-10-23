#include "RSTfile.hpp"

RSTfile::RSTfile(std::string filename,std::string path):
	RST(),
	path_(path),
	filename_(filename),
	create_pdf_(false)
{ }

RSTfile::~RSTfile() { 
	Write w(path_ + filename_ + ".rst");
	rst += RST_np;
	w<<rst;
	Linux command;
	command("rst2html " + path_  + filename_ + ".rst " + path_ + filename_ + ".html");  
	if(command.status()){
		std::cerr<<"RSTfile::~RSTfile() : the command function called returns an error for the html creation"<<std::endl; 
	} else {
		if(create_pdf_){
			command("rst2latex " + path_ + filename_ + ".rst " + path_ + filename_ + ".tex"); 
			command("pdflatex -output-directory " + path_ + " "  + filename_ + ".tex");
			if(command.status()){std::cerr<<"RSTfile::~RSTfile() : the command function called returns an error for the pdf creation"<<std::endl; }
		}
	}
}

