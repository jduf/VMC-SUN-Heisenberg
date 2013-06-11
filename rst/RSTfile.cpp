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
	std::string command("rst2html " + path  + filename + ".rst " + path + filename + ".html");  
	int status(0);
	status = system(command.c_str());
	if(status != 0){
		std::cerr<<"RSTfile::~RSTfile() : the command function called returns an error for the html creation"<<std::endl; 
	} else {
		if(create_pdf){
			command="rst2latex " + path + filename + ".rst " + path + filename + ".tex";  
			status = system(command.c_str());
			command = "pdflatex -output-directory " + path + " "  + filename + ".tex"; 
			status = system(command.c_str());
			if(status != 0){std::cerr<<"RSTfile::~RSTfile() : the command function called returns an error for the pdf creation"<<std::endl; }
		}
	}
}

