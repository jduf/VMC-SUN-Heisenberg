#include "RSTfile.hpp"

RSTfile::RSTfile(std::string filename,std::string path):
	RST(),
	path(path),
	filename(filename),
	w(path + filename + ".rst"),
	create_pdf(false)
{ }

RSTfile::~RSTfile() { 
	rst += RST_np;
	w<<rst;
	Linux command;
	command("rst2html " + path  + filename + ".rst " + path + filename + ".html");  
	if(command.status()){
		std::cerr<<"RSTfile::~RSTfile() : the command function called returns an error for the html creation"<<std::endl; 
	} else {
		if(create_pdf){
			command("rst2latex " + path + filename + ".rst " + path + filename + ".tex"); 
			command("pdflatex -output-directory " + path + " "  + filename + ".tex");
			if(command.status()){std::cerr<<"RSTfile::~RSTfile() : the command function called returns an error for the pdf creation"<<std::endl; }
		}
	}
}

