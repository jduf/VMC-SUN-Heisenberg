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
	Vector<double> vec;
	vec = P.get<Vector<double> >("vec");
	std::cout<<vec<<std::endl;
	if(P.is_vector("param")){
		Vector<double> param(P.get<Vector<double> >("param"));
		std::cout<<param<<std::endl;
	} else {
		double param(P.get<double>("param"));
		std::cout<<param<<std::endl;
	}
}
