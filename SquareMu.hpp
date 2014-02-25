#ifndef DEF_SQUAREMU
#define DEF_SQUAREMU

#include "Square.hpp"

class SquareMu: public Square<double>{
	public:
		SquareMu(Container const& param);
		~SquareMu();

		void create(double mu);
		void study();
		void save();
		void get_param(Container& param);

	protected:
		double mu_;

		void compute_T(unsigned int alpha);
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif

