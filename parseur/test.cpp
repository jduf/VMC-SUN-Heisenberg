#include "Parseur.hpp"

#include <string>

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	double a(0);
	std::string b("");
	double c(0);
	unsigned int d(0);
	P.set("a",a);
	P.set("b",b);
	P.set("1",c);
	d = P.get<unsigned int>("0");
	std::cout<<"a="<<a<<std::endl;
	std::cout<<"b="<<b<<std::endl;
	std::cout<<"c="<<c<<std::endl;
	std::cout<<"d="<<d<<std::endl;
	
}



