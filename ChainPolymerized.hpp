#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(Container const& param);
		~ChainPolymerized();

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
