#include "Analyse.hpp"

int main(int argc, char* argv[]){
	if(argc==1){
		Analyse a("");
	} else {
		Analyse a(argv[1]);
	}
}
