#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Matrice.hpp"

class System{
	public:
		System(unsigned int N_spin, unsigned int N_m, unsigned int dim);
		~System();

		unsigned int const N_spin, N_m, N_site;
		unsigned int Nx, Ny, dim;
		unsigned int *nts;
		Matrice U;

	private:
		System();
		void create_U(unsigned int dim);
		void create_nts(unsigned int dim);
};
#endif
