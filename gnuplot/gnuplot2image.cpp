/*! @file gnuplot2image.cpp*/

#include "Linux.hpp"
#include "Miscellaneous.hpp"

int main(int argc,char* argv[]){
	if(argc!=1){
		std::string filename(argv[1]);
		const std::string ext(".gp");
		if ( filename != ext && filename.size() > ext.size() && filename.substr(filename.size() - ext.size()) == ".gp" ) {
			filename = filename.substr(0, filename.size() - ext.size());
			std::string path("");

			std::vector<std::string> tmp(my::string_split(path+filename,'/'));
			for(unsigned int i(0);i<tmp.size()-1;i++){ path += tmp[i] + "/"; }
			filename = tmp[tmp.size()-1];

			Linux command;
			if(path[0] != '/'){ path = command.pwd() + path; }
			std::string texfile("gnuplot2image_tmp");
			command(Linux::gp2latex("/tmp/"+texfile,path,filename));
			if(!command.status()){
				command(Linux::pdflatex("/tmp/",texfile),true);
				command(Linux::pdf2png("/tmp/" + texfile, path + filename));
				command("mv /tmp/" + texfile + ".pdf " + path + filename + ".pdf");
				command("rm /tmp/" + texfile + "*");
			} else {
				std::cerr<<"gnuplot : the script contains mistakes"<<std::endl;
			}
		} else { std::cerr<<"gnuplot2pdf : the filename must have a '.gp' extension"<<std::endl; }
	} else { std::cerr<<"gnuplot2pdf : a gnuplot file (*.gp) must be given as argument"<<std::endl; }
}
