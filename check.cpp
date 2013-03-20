#include "Read.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"

#include <string>
void check(std::string filename);

int main(int argc, char* argv[]){
	if(argc!=2){
		std::cerr<<"check : wrong number of arguements"<<std::endl;
	} else {
		check(argv[1]);
	}
}

void check(std::string filename){
	std::cout<<filename<<std::endl;
	Read r(filename.c_str());
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
