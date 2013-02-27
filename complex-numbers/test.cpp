#include "Complex.hpp"
using namespace std;

int main(){
	Complex c1;
	Complex c2(2,1);
	Complex c3(3,5);
	Complex c4(3,6);
	std::cout<<exp(ii(M_PI)).ProjRe()<<std::endl;
	std::cout<<exp(ii(M_PI/4.0)).ProjRe()<<std::endl;
	std::cout<<exp(ii(M_PI/2.0))<<std::endl;
	std::cout<<exp(ii(3.0*M_PI/2.0))<<std::endl;
	std::cout<<exp(ii(2.0*M_PI))<<std::endl;
	std::cout<<exp(c4)<<std::endl;


	std::cout<<c2<<std::endl;
	std::cout<<c2/2.0<<std::endl;
	std::cout<<2.0/c2<<std::endl;
	std::cout<<c3/c2<<std::endl;
}
