/*!  @file jdread.cpp */
#include "IOFiles.hpp"

int main(int argc, char* argv[]){
	if(argc == 2){
		IOFiles r(argv[1],false,false);
		std::cout<<r.get_header()<<std::endl;
	} else {
		std::cerr<<"jdread : takes exactly one input argument : ./jdread filename"<<std::endl;
	}
}	
