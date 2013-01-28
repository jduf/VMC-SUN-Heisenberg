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
		State(System *S);
		//cérateurs de copie profonde (néessaire car s et wis sont des pointeurs)
		State(unsigned int N_m, unsigned int N_spin);
		State(State const& s);
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
		double* get_mat(unsigned int i) const;
		unsigned int get_changed_column(unsigned int i) const;

	private:
		// pas de créateur par défaut, car pas d'allocation mémoire
		State();

		System *S;
		Matrice *A;
		Matrice *Ainv;
		unsigned int *s;
		unsigned int *wis;
		unsigned int column_changed[2];
		unsigned int spin_changed[2];
		double det;

		void init_A(unsigned int N_m, unsigned int N_spin);
		void compute_matrices();
		void compute_det(); // sign(permutation) => ratio det tjrs - ??? 
		void modify_matrix(unsigned int a, unsigned int b);
};

std::ostream& operator<<(std::ostream& flux, State const& S);
double operator/(State const& Snew, State const& Sold);
double compute_ratio_det(Matrice const& mold, Matrice const& mnew, unsigned int col, unsigned int N);
#endif
