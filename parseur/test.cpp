#include "Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	double a(0);
	std::string b("");
	double c(0);
	unsigned int d(0);
	P.get("a",a);
	P.get("b",b);
	P.get("1",c);
	d = P.get<unsigned int>("0");
	std::cout<<"a="<<a<<std::endl;
	std::cout<<"b="<<b<<std::endl;
	std::cout<<"c="<<c<<std::endl;
	std::cout<<"d="<<d<<std::endl;
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
