#include "CreateSystem.hpp"
#include "Parseur.hpp"

#include <sstream>
#include <string>

void check(std::string lattice);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex(false);
	std::string lattice;

	P.set("N_spin",N_spin);	
	P.set("N_m",N_m);	
	P.set("N_n",N_n);
	P.set("complex",is_complex);
	switch(N_n){
		case 2:
			{
				lattice="chain";
				break;
			}
		case 3:
			{
				lattice="honeycomb";
				break;
			}
		case 4:
			{
				lattice="square";
				break;
			}
		default:
			{
				std::cerr<<"main : lattice type undefined"<<std::endl;
			}
	}
	std::stringstream ss1;
	std::stringstream ss2;
	ss1<<N_spin;
	ss2<<N_spin*N_m;
	lattice += "-N"+ss1.str() + "-S"+ss2.str();
	if(is_complex){ lattice += "-1"; }
	else{ lattice += "-0"; }

	if(is_complex){
		CreateState<std::complex<double> > CS(N_m,N_spin,N_n,lattice);
	} else {
		CreateState<double> CS(N_m,N_spin,N_n,lattice);
	}

	check(lattice);
}

void check(std::string lattice){
	std::cout<<lattice<<std::endl;
	Read r(lattice.c_str());
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex;

	r>>is_complex>>N_spin>>N_m>>N_n;
	Array2D<unsigned int> nts(N_spin*N_n*N_m/2,2);
	r>>nts;
	std::cout<<"N_spin="<<N_spin<<" N_m="<<N_m<<" N_n="<<N_n<<std::endl;
	std::cout<<"nts="<<std::endl;
	std::cout<<nts<<std::endl;
	if(is_complex){
		Matrice<std::complex<double> > H;
		Matrice<std::complex<double> > EVec;
		r>>H>>EVec;
		std::cout<<"H="<<std::endl;
		std::cout<<H<<std::endl;
		std::cout<<"EVec="<<std::endl;
		std::cout<<EVec<<std::endl;
	} else {
		Matrice<double> H;
		Matrice<double> EVec;
		r>>H>>EVec;
		std::cout<<"H="<<std::endl;
		std::cout<<H<<std::endl;
		std::cout<<"EVec="<<std::endl;
		std::cout<<EVec<<std::endl;
	}
}
