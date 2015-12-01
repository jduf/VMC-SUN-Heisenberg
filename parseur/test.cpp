#include "Vector.hpp"
#include "Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);

	double a(0);
	std::string b("");
	unsigned int c(0);
	P.get("a",a);
	P.get("b",b);
	c = P.get<unsigned int>("c");
	std::cout<<"a="<<a<<std::endl;
	std::cout<<"b="<<b<<std::endl;
	std::cout<<"c="<<c<<std::endl;
	std::cout<<Vector<double>(P.get<std::vector<double> >("s_vec"))<<std::endl;
	std::cout<<Vector<double>(P.get<std::vector<double> >("vec"))<<std::endl;
}
