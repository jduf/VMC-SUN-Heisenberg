#include "Read.hpp"
#include "Write.hpp"
#include "Array2D.hpp"
#include "Lapack.hpp"
#include "Matrice.hpp"

#include <complex>

template<typename M>
class CreateState{
	public:
		CreateState(unsigned int N_m, unsigned int N_spin, unsigned int N_n, std::string filename);

	private:
		std::string const filename;
		char mat_type;
		Write w;
		unsigned int const N_m, N_spin, N_n, N_site;
		Matrice<M> H;

		void init();
		void compute_EVec();

		void init_sts();
};

template<typename M>
CreateState<M>::CreateState(unsigned int N_m, unsigned int N_spin, unsigned int N_n, std::string filename):
	filename(filename),
	mat_type('U'),
	w(),
	N_m(N_m),
	N_spin(N_spin), 
	N_n(N_n),
	N_site(N_m*N_spin),
	H(N_spin*N_m)
{
	init();
	init_sts();
	compute_EVec();
}


template<typename M>
void CreateState<M>::init_sts(){
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
	w<<nts<<H;
}

template<typename M>
void CreateState<M>::compute_EVec(){
	Lapack<M> ES(H.ptr(),H.size(), mat_type);
	Vecteur<double> EVal(H.size());
	ES.eigensystem(EVal);
	w<<H;
}
