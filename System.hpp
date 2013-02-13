#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Vecteur.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"
#include "Read.hpp"
#include <cstdlib>

class System{
	public:
		System(unsigned int N_spin, unsigned int N_m, unsigned int N_n, std::string filename);
		~System();

		unsigned int const N_spin, N_m, N_n, N_site;
		unsigned int *nts;
		
		void update_state();
		void swap();
		void swap(unsigned int a, unsigned int b);
		double compute_ratio();
		double det();

		unsigned int const& operator[](unsigned int const& i) const { return wis[i];};

		void print();

	private:
		System();
		Matrice<double> *A;
		Matrice<double> *Ainv;
		Matrice<double> tmp_mat;
		unsigned int *s;
		unsigned int *wis;
		unsigned int cc[2]; // column changed
		unsigned int mc[2]; // matrix changed
		double w[2];

		void create_nts(Matrice<double> const& U);
		void init_state(Matrice<double> const& U);
};

std::ostream& operator<<(std::ostream& flux, System const& S);
#endif
