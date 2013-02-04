#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Matrice.hpp"
#include "Lapack.hpp"
#include <cstdlib>

class System{
	public:
		System(unsigned int N_spin, unsigned int N_m, unsigned int dim);
		~System();

		unsigned int const N_spin, N_m, N_site, dim;
		unsigned int *nts;
		
		void update_state();
		void swap();
		void swap(unsigned int a, unsigned int b);
		double compute_ratio();
		double det();

		unsigned int operator[](unsigned int i) const;

		void print();

	private:
		System();
		Matrice *A;
		Matrice *Ainv;
		unsigned int *s;
		unsigned int *wis;
		unsigned int Nx,Ny;
		unsigned int cc[2]; // column changed
		unsigned int mc[2]; // matrix changed
		double w[2];

		void create_U(Matrice& U, unsigned int dim);
		void create_nts(unsigned int dim);
		void init_state(Matrice const& U);
};
#endif
