#include "Gnuplot.hpp"
#include "Parseur.hpp"

#include<iostream>

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	std::string s(P.get<std::string>("0"));
	if(!P.status()){
		const std::string ext(".gp");
		if ( s != ext && s.size() > ext.size() && s.substr(s.size() - ext.size()) == ".gp" ) {
			s = s.substr(0, s.size() - ext.size());
			std::string path;
			if(s[0] != '/'){
				Linux c;
				path = c.pwd();
			}

			std::vector<std::string> tmp(string_split(path+s,'/'));
			path = "";
			for(unsigned int i(0);i<tmp.size()-1;i++){ path += tmp[i] + "/"; }
			s = tmp[tmp.size()-1];

			std::cout<<path<<" "<<s<<std::endl;

			Gnuplot gp(path,s);
			gp.create_image(false);
		} else { std::cerr<<"gnuplot2pdf : the filename must have a '.gp' extension"<<std::endl; }
	}
}
