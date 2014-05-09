#include "RSTFile.hpp"

RSTFile::RSTFile(std::string path,std::string filename):
	RST(),
	path_(path),
	filename_(filename),
	pdf_(false)
{}

RSTFile::~RSTFile() { 
	IOFiles w(path_ + filename_ + ".rst",true);
	rst += RST_np;
	w<<rst;
	Linux command;
	//command("rst2html --stylesheet=/home/jdufour/travail/cpp-dev/rst/css/voidspace.css --field-name-limit=0 " + path_  + filename_ + ".rst " + path_ + filename_ + ".html");  
	command("rst2html --field-name-limit=0 " + path_  + filename_ + ".rst " + path_ + filename_ + ".html");  
	if(command.status()){
		std::cerr<<"RSTFile::~RSTFile() : the command function called returns an error for the html creation"<<std::endl; 
	} else {
		if(pdf_){
			command("rst2latex " + path_ + filename_ + ".rst " + path_ + filename_ + ".tex"); 
			command("pdflatex -output-directory " + path_ + " "  + filename_ + ".tex");
			if(command.status()){std::cerr<<"RSTFile::~RSTFile() : the command function called returns an error for the pdf creation"<<std::endl; }
		}
	}
}

