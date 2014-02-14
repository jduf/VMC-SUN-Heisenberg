#ifndef DEF_CHAINTRIMERIZED
#define DEF_CHAINTRIMERIZED

#include "Chain.hpp"

class ChainTrimerized: public Chain<double> {
	public:
		ChainTrimerized(Container const& param);
		~ChainTrimerized();

		void save();
		void create(Vector<double> const& t);
		void properties(Container& c);

	private:
		void compute_P(Matrix<double>& P);
		void compute_T(Vector<double> const& t);

		Vector<double> t_;
};
#endif
