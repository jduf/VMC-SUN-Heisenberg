#ifndef DEF_STATE
#define DEF_STATE

#include<vector>
#include<cstdlib>
#include<iostream>
#include<armadillo>

class State{
	public:
		//créateur initial prenant un pointeur vers une arma::mat
		State(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U);
		//cérateurs de copie profonde (néessaire car s et wis sont des pointeurs)
		State(State const& s);
		State(unsigned int N_site);
		//destructeur necessaire car utilisation de new dans certains créateurs
		~State();
		//affectation profonde (pointeurs)
		State& operator=(State const& s);

		// échange deux site de spin différent
		State swap() const;
		// échange les site a et b (pouvant être de même spin)
		State swap(unsigned int a, unsigned int b) const;
		void print() const;
		inline double Det() const {return det;}

	private:
		State();
		unsigned int N_spin, N_m;
		unsigned int *s;
		unsigned int *wis;
		double det;
		std::vector<arma::Mat<double> > A;
		arma::Mat<double> *U; 

		void compute_matrices();
		void compute_det(); // sign(permutation) => ratio det tjrs - ??? 
};
#endif
