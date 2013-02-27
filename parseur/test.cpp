#include "Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	unsigned int N(0);
	P.set("N",N);
	std::cout<<N<<std::endl;
	double a(0);
	P.set("a",a);
	std::cout<<a<<std::endl;
}



