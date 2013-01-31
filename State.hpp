#ifndef DEF_STATE
#define DEF_STATE

#include <cstdlib>
#include <iostream>
#include "System.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"

class State{
	public:
		//créateur initial prenant un pointeur vers une arma::mat
		State(System *S,bool whole_copy);
		//destructeur necessaire car utilisation de new dans certains créateurs
		~State();
		//affectation profonde (pointeurs)
		State& operator=(State const& s);

		// échange deux site de spin différent
		State swap() const;
		// échange les site a et b (pouvant être de même spin)
		State swap(unsigned int a, unsigned int b) const;
		void print() const;
		void color(std::ostream& flux) const;
		double divided() const ;
		double Det() { return det;};

	private:
		// pas de créateur par défaut, car pas d'allocation mémoire
		State();
		State(State const& s); //copie de surface

		System *S;
		unsigned int *s;
		unsigned int *wis;
		unsigned int cc[2]; // column changed
		unsigned int mc[2]; // matrix changed
		double det;
		bool delete_matrix;
		bool whole_copy;

		void compute_det(); // sign(permutation) => ratio det tjrs - ??? 
		void init_matrices(unsigned int N_m, unsigned int N_spin);
};

std::ostream& operator<<(std::ostream& flux, State const& S);
double operator/(State const& Snew, State const& Sold);
#endif
