#include "Read.hpp"
#include "Write.hpp"
#include "Array2D.hpp"
#include "Lapack.hpp"
#include "Matrice.hpp"

#include <complex>
#include <sstream>

template<typename M>
class CreateState{
	public:
		CreateState(unsigned int N_m, unsigned int N_spin, unsigned int N_n, std::string filename);
		~CreateState();

	private:
		std::string filename;
		std::string filename_comp;
		bool is_complex;
		char mat_type;
		Write w;
		unsigned int const N_m, N_spin, N_n, N_site;
		Matrice<M> H,EVec;
		Array2D<unsigned int> sts;

		void compute_H();
		void compute_EVec();
		void compute_sts();
};

template<typename M>
CreateState<M>::CreateState(unsigned int N_m, unsigned int N_spin, unsigned int N_n, std::string filename):
	filename(filename),
	filename_comp(""),
	is_complex(false),
	mat_type('U'),
	w(),
	N_m(N_m),
	N_spin(N_spin), 
	N_n(N_n),
	N_site(N_m*N_spin),
	H(N_spin*N_m),
	EVec(N_spin*N_m),
	sts(N_spin*N_m*N_n/2,2)
{
	compute_H();
	compute_sts();
	compute_EVec();
}


template<typename M>
CreateState<M>::~CreateState(){
	w.open(filename+filename_comp);
	w<<is_complex<<N_spin<<N_m<<N_n<<sts<<H<<EVec;
}

template<typename M>
void CreateState<M>::compute_sts(){
	unsigned int k(0);
	for(unsigned int i(0); i<N_site;i++){
		for(unsigned int j(i+1); j<N_site;j++){
			if ( std::abs(H(i,j)) > 1e-4){
				sts(k,0) = i;
				sts(k,1) = j;
				k++;
			}
		}
	}
}

template<typename M>
void CreateState<M>::compute_EVec(){
	EVec = H;
	Lapack<M> ES(EVec.ptr(),N_site, mat_type);
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
	if(std::abs(EVal(N_m) - EVal(N_m-1))<1e-6){
		filename_comp = "-VPDNF";
	}
}
