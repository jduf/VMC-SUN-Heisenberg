#ifndef DEF_CHAINSAS
#define DEF_CHAINSAS

#include "Chain.hpp"

/*!{For k=2, enforce a symmetric irrep by tuning J/J'*/
class ChainSAS: public Chain<double> {
	public:
		ChainSAS(System const& s, Vector<double> const& t);
		~ChainSAS() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_;//!< polymerization parameter

		void compute_H(unsigned int const& c);
		void display_results(){};

		unsigned int set_spuc(unsigned int const N);
};
#endif
