#ifndef DEF_TRIANGLEMU
#define DEF_TRIANGLEMU

#include "Triangle.hpp"

class TriangleMu: public Triangle<double>{
	public:
		TriangleMu(Container const& param);
		~TriangleMu();

		void create(double mu);
		void study();
		void save();
		void get_param(Container& param);

	protected:
		double mu_;

		void compute_T(unsigned int alpha);
		void compute_P();
		void lattice();
};
#endif

