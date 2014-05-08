/*! @file gnuplot2image.cpp*/

#include "Gnuplot.hpp"
#include "Parseur.hpp"

#include<iostream>

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	std::string filename(P.get<std::string>("0"));
	if(!P.status()){
		const std::string ext(".gp");
		if ( filename != ext && filename.size() > ext.size() && filename.substr(filename.size() - ext.size()) == ".gp" ) {
			filename = filename.substr(0, filename.size() - ext.size());
			std::string path;
			if(filename[0] != '/'){
				Linux c;
				path = c.pwd();
			}

			std::vector<std::string> tmp(string_split(path+filename,'/'));
			path = "";
			for(unsigned int i(0);i<tmp.size()-1;i++){ path += tmp[i] + "/"; }
			filename = tmp[tmp.size()-1];

			Gnuplot gp(path,filename);
			gp.create_image(true);
		} else { std::cerr<<"gnuplot2pdf : the filename must have a '.gp' extension"<<std::endl; }
	}
}
