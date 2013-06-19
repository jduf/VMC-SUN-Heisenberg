#include "Parseur.hpp"

#include <string>

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	double a(0);
	std::string b("");
	double c(0);
	unsigned int d(0);
	P.set("a",a);
	P.print();
	std::cout<<std::endl;
	P.set("b",b);
	P.print();
	std::cout<<std::endl;
	P.set("1",c);
	P.print();
	std::cout<<std::endl;
	d = P.get<unsigned int>("0");
	P.print();
	std::cout<<std::endl;
	std::cout<<"a="<<a<<std::endl;
	std::cout<<"b="<<b<<std::endl;
	std::cout<<"c="<<c<<std::endl;
	std::cout<<"d="<<d<<std::endl;
	
}



