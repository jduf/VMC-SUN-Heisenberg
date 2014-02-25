#ifndef DEF_CHAINDIMERIZED
#define DEF_CHAINDIMERIZED

#include "Chain.hpp"

class ChainDimerized: public Chain<double> {
	public:
		ChainDimerized(Container const& param);
		~ChainDimerized();

		void create(double delta);
		void study();
		void save();
		void get_param(Container& param);

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		double delta_;
};
#endif
