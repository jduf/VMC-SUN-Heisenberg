/*! @file gnuplot2image.cpp*/

#include "Linux.hpp"
#include "Miscellaneous.hpp"

int main(int argc,char* argv[]){
	const std::string ext(".gp");
	if(argc==2 || argc==3){
		std::string filename(argv[1]);
		if ( filename != ext && filename.size() > ext.size() && filename.substr(filename.size() - ext.size()) == ext ) {
			filename = filename.substr(0, filename.size() - ext.size());
			std::string path("");

			std::vector<std::string> tmp(my::string_split(path+filename,'/'));
			for(unsigned int i(0);i<tmp.size()-1;i++){ path += tmp[i] + "/"; }
			filename = tmp[tmp.size()-1];

			Linux command;
			if(path[0] != '/'){ path = command.pwd() + path; }
			std::string texfile("/tmp/gp2latex_tmp");
			command(Linux::gp2latex(texfile,path,filename),false);
			if(!command.status()){
				command(Linux::pdflatex("/tmp/",texfile),true);
				if(argc==3){
					std::string option(argv[2]);
					if(option=="png"){ command(Linux::pdf2png(texfile, path + filename),true); }
					else { std::cerr<<__PRETTY_FUNCTION__<<" : unknown option (only possible is 'png')"<<std::endl; }
				}
				command("mv " + texfile + ".pdf " + path + filename + ".pdf",false);
				command("rm " + texfile + "*",false);
			} else { std::cerr<<__PRETTY_FUNCTION__<<" : the script contains mistakes"<<std::endl; }
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : the filename must have a '"<<ext<<"' extension"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : a Gnuplot file (*"<<ext<<") must be given as argument (will generate a 'png' file if the second argument is 'png') "<<std::endl; }
}
