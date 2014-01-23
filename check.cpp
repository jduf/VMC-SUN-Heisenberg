/*!  @file check.cpp */

#include "ExtractSystem.hpp"

int main(int argc, char* argv[]){
	if(argc!=2){
		std::cerr<<"check : wrong number of arguements"<<std::endl;
	} else {
		ExtractSystem system(argv[1]);
		system.extract();
		system.print();
	}
}
