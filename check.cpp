#include "Read.hpp"
#include "Matrix.hpp"

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
	std::cout<<r.get_header()<<std::endl;
	unsigned int N_spin(0),N_row(0),N_col(0),N_m(0);
	double bc(0.0);
	bool is_complex;
	Matrix<unsigned int> sts;

	r>>is_complex>>N_spin>>N_m>>sts;
	for(unsigned int i(0);i<sts.row();i++){
		std::cout<<sts(i,0)<<" "<<sts(i,1)<<std::endl;
	}
	if(is_complex){
		Matrix<std::complex<double> > EVec;
		r>>EVec;
		EVec.chop();
		std::cout<<"EVec="<<std::endl;
		std::cout<<EVec<<std::endl;
	} else {
		Matrix<double> EVec;
		r>>EVec;
		EVec.chop();
		std::cout<<"EVec="<<std::endl;
		std::cout<<EVec<<std::endl;
	}
	r>>bc>>N_row>>N_col;	
	std::cout<<"N_spin="<<N_spin
		<<" N_site="<<N_m*N_spin
		<<" N_row="<<N_row
		<<" N_col="<<N_col
		<<" bc="<<bc<<std::endl;
}
