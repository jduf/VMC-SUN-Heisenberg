#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Matrice.hpp"
#include "Lapack.hpp"

class System{
	public:
		System(unsigned int N_spin, unsigned int N_m, unsigned int dim);
		~System();

		unsigned int const N_spin, N_m, N_site, dim;
		unsigned int *nts;
		Matrice U;
		Matrice *A;
		Matrice *Ainv;
		double w[2];

		void update_matrices(unsigned int const mc[], unsigned int const cc[]);

	private:
		System();
		unsigned int Nx,Ny;
		void create_U(unsigned int dim);
		void create_nts(unsigned int dim);
};
#endif
