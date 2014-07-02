/*!  @file jdread.cpp */
#include "IOFiles.hpp"
#include "Parseur.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string filename(P.get<std::string>("0"));
	if(!P.status()){
		IOFiles r(filename,false);
		std::cout<<r.get_header()<<std::endl;
	} else {
		std::cerr<<"jdread : take exactly one input argument"<<std::endl;
		std::cerr<<"usage -> ./jdread filename"<<std::endl;
		std::cerr<<"note  -> the file should be given without the .jdbin extension"<<std::endl;
	}
}	
