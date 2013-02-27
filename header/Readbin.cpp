#include "Read.hpp"
#include "Header.hpp"

#include <string>
#include <iostream>

int main(int argc, char* argv[]){
	if(argc==2){
		std::string filename(argv[1]);
		Read r(filename);
		Header h(r);
		h.show();
	} else {
		std::cerr<<"readbin : take exactly one input argument"<<std::endl;
		std::cerr<<"usage -> ./readbin filename"<<std::endl;
		std::cerr<<"note  -> the file should be given without the .bin extension"<<std::endl;
	}
}	
