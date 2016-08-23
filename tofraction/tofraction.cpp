#include "Miscellaneous.hpp"

int main(int argc, char* argv[]){
	if(argc==2){
		double d(my::string2type<double>(argv[1]));
		unsigned long long num(1);
		unsigned long long den(1);
		double sign;
		unsigned int iter(my::to_fraction(d,num,den,sign));
		std::cout<<"result obtained in "<<iter<<" iterations : "<<d<<" = "<<sign*num/den<<" = "<<sign*num<<"/"<<den<<std::endl;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : to first argument must be a real scalar, other argument are ignored"<<std::endl; }
}
