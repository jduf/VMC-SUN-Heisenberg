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
	unsigned int N(0),m(0);
	double bc(0.0);
	std::string wf;
	Matrix<unsigned int> sts;
	r>>wf>>N>>m>>sts;
	std::cout<<"sts="<<std::endl;
	for(unsigned int i(0);i<sts.row();i++){
		std::cout<<sts(i,0)<<" "<<sts(i,1)<<std::endl;
	}
	if( wf == "chain" ){
		Matrix<double> EVec;
		r>>EVec>>bc;
		std::cout<<"N_spin="<<N
			<<" N_site="<<m*N
			<<" bc="<<bc
			<<std::endl;
	}
	if( wf == "fermi" ){
		unsigned int Lx(0),Ly(0);
		Matrix<double> EVec;
		r>>EVec>>bc>>Lx>>Ly;
		std::cout<<"N_spin="<<N
			<<" N_site="<<m*N
			<<" Lx="<<Lx
			<<" Ly="<<Ly
			<<" bc="<<bc
			<<std::endl;
	}
	if( wf == "mu" ){
		unsigned int Lx(0),Ly(0);
		double mu(0);
		Matrix<double> EVec;
		r>>EVec>>bc>>Lx>>Ly>>mu;
		std::cout<<"N_spin="<<N
			<<" N_site="<<m*N
			<<" Lx="<<Lx
			<<" Ly="<<Ly
			<<" bc="<<bc
			<<" mu="<<mu
			<<std::endl;
		std::cout<<EVec.chop()<<std::endl;
	}
	if( wf == "csl" ){
		unsigned int Lx(0),Ly(0);
		Matrix<std::complex<double> > EVec;
		r>>EVec>>bc>>Lx>>Ly;
		std::cout<<"N_spin="<<N
			<<" N_site="<<m*N
			<<" Lx="<<Lx
			<<" Ly="<<Ly
			<<" bc="<<bc
			<<std::endl;
	}
	if( wf == "honeycomb" ){
		unsigned int Lx(0),Ly(0);
		Matrix<double> EVec;
		r>>EVec>>bc>>Lx>>Ly;
		std::cout<<"N_spin="<<N
			<<" N_site="<<m*N
			<<" Lx="<<Lx
			<<" Ly="<<Ly
			<<" bc="<<bc
			<<std::endl
			<<" U="
			<<std::endl
			<<EVec
			<<std::endl;
	}
}
