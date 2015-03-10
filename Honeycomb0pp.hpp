#ifndef DEF_HONEYCOMB0PP
#define DEF_HONEYCOMB0PP

#include "Honeycomb.hpp"

class Honeycomb0pp: public Honeycomb<double>{
	public:
		Honeycomb0pp(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, double td);
		~Honeycomb0pp(){}

		void create();
		void save() const;
		void check();

	protected:
		double td_;

		void compute_H();
		void lattice();
		std::string extract_level_7();
		std::string extract_level_6();

		Matrix<double> set_ab(){
			Matrix<double> tmp(2,2);
			tmp(0,0) = 1;
			tmp(1,0) = 1;
			tmp(0,1) = 2;
			tmp(1,1) = -1;
			return tmp;
		}
};
#endif

