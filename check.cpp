#include "Read.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"

#include <string>
#include <iostream>
void check(std::string filename);

int main(int argc, char* argv[]){
	if(argc!=2){
		std::cerr<<"check : wrong number of arguements"<<std::endl;
	} else {
		check(argv[1]);
	}
}

void check(std::string filename){
	Read r(filename);
	std::cout<<r.header()<<std::endl;
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex;

	r>>is_complex>>N_spin>>N_m>>N_n;
	Array2D<unsigned int> sts(N_spin*N_n*N_m/2,2);
	Matrice<double> H(N_m*N_spin);
	r>>sts>>H;
	std::cout<<"N_spin="<<N_spin<<" N_m="<<N_m<<" N_n="<<N_n<<std::endl;
	std::cout<<"nts="<<std::endl;
	std::cout<<sts<<std::endl;
	std::cout<<"H="<<std::endl;
	std::cout<<H<<std::endl;
	if(is_complex){
		Matrice<std::complex<double> > EVec(N_m*N_spin);
		r>>EVec;
		std::cout<<"EVec="<<std::endl;
		std::cout<<EVec<<std::endl;
	} else {
		Matrice<double> EVec(N_m*N_spin);
		r>>EVec;
		std::cout<<"EVec="<<std::endl;
		std::cout<<EVec<<std::endl;
	}
}
