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
	command("rst2html --stylesheet=/home/jdufour/travail/cpp-dev/rst/css/nice.css --field-name-limit=0 " + path_  + filename_ + ".rst " + path_ + filename_ + ".html");  
	//command("rst2html --stylesheet=/home/jdufour/rst2html5.css --field-name-limit=0 " + path_  + filename_ + ".rst " + path_ + filename_ + ".html");  
	//command("rst2html --field-name-limit=0 " + path_  + filename_ + ".rst " + path_ + filename_ + ".html");  
	if(command.status()){
		std::cerr<<"RSTFile::~RSTFile() : the command function called returns an error for the html creation"<<std::endl; 
	} else {
		if(pdf){
			command("rst2latex " + path_ + filename_ + ".rst " + path_ + filename_ + ".tex"); 
			command("pdflatex -output-directory " + path_ + " "  + filename_ + ".tex");
			if(command.status()){std::cerr<<"RSTFile::~RSTFile() : the command function called returns an error for the pdf creation"<<std::endl; }
		}
	}
}

