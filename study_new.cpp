
#include "Analyse.hpp"

int main(int argc, char* argv[]){
	
	std::cout<<argc<<std::endl;
	if(argc==1){
		Analyse a("");
	} else {
		Analyse a(argv[1]);
	}
}
