#include "Write.hpp"
#include "Read.hpp"
#include "Parseur.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"

void compute_EVec(Matrice<double>& H, Write& w);

void init_sts(Matrice<double> const& H, unsigned int N_n, Write& w);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string sysname("hopmat24");;
	unsigned int N_spin(4), N_m(6), N_n(3);

	Read r(sysname.c_str(),false);
	Matrice<double> H(24);
	r>>H;
	Write w(sysname.c_str());
	bool is_complex(false);
	w<<is_complex<<N_spin<<N_m<<N_n;
	init_sts(H,N_n,w);
	compute_EVec(H,w);
}

void compute_EVec(Matrice<double>& H, Write& w){
	Lapack<double> ES(H.ptr(),H.size(),'S');
	Vecteur<double> EVal(H.size());
	ES.eigensystem(EVal);
	EVal.print();
	H.print();
	w<<H;
}

void init_sts(Matrice<double> const& H, unsigned int N_n, Write& w){
	unsigned int k(0),N_site(H.size());
	Array2D<unsigned int> nts(N_site*N_n/2,2);
	for(unsigned int i(0); i<N_site;i++){
		for(unsigned int j(i+1); j<N_site;j++){
			if ( std::abs(H(i,j)) > 1e-4){
				nts(k,0) = i;
				nts(k,1) = j;
				k++;
			}
		}
	}
	w<<nts;
}
