#ifndef DEF_CHAINTRIMERIZED
#define DEF_CHAINTRIMERIZED

#include "Chain.hpp"

class ChainTrimerized: public Chain<double> {
	public:
		ChainTrimerized(Container const& param);
		~ChainTrimerized();

		void create(Vector<double> const& t);
		void study();
		void save();
		void get_param(Container& param);

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		Vector<double> t_;
};
#endif
