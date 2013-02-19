#include "Read.hpp"
#include "Write.hpp"
#include "Parseur.hpp"
#include "Array2D.hpp"
#include "Lapack.hpp"
#include "Matrice.hpp"

#include <sstream>
#include <complex>

void init_H(Matrice<double>& H, unsigned int N_n);
void init_H(Matrice<std::complex<double> >& H, unsigned int N_n);
void compute_EVec(Matrice<double>& H, Write& w);
void compute_EVec(Matrice<std::complex<double> >& H, Write& w);

template<typename M>
void init_sts(Matrice<M> const& H, unsigned int N_n, Write& w);

void check(std::string lattice);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex(false);
	std::string lattice;

	P.set("N_spin",N_spin);	
	P.set("N_m",N_m);	
	P.set("lattice",lattice);	
	P.set("complex",is_complex);
	if(lattice=="chain"){ N_n = 2;}
	if(lattice=="honeycomb"){ N_n = 3;}
	if(lattice=="square"){ N_n = 4;}
	std::stringstream ss1;
	std::stringstream ss2;
	ss1<<N_spin;
	ss2<<N_spin*N_m;
	lattice += "-N"+ss1.str() + "-S"+ss2.str();
	if(is_complex){ lattice += "-1"; }
	else{ lattice += "-0"; }
	//N_spin = 3; N_m = 4; N_n = 4; lattice="check"; is_complex=false;
	std::cout<<lattice<<std::endl;

	Write w(lattice.c_str());
	w<<is_complex<<N_spin<<N_m<<N_n;
	std::cout<<is_complex<<std::endl;
	if(is_complex){
		Matrice<std::complex<double> > H(N_spin*N_m);
		init_H(H,N_n);
		std::cout<<"Hcomp"<<std::endl;
		std::cout<<H<<std::endl;
		init_sts(H,N_n,w);
		compute_EVec(H,w);
	} else {
		Matrice<double> H(N_spin*N_m);
		init_H(H,N_n);
		init_sts(H,N_n,w);
		compute_EVec(H,w);
	}

	check(lattice);
}
