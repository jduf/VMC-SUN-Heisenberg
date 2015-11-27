#ifndef DEF_TRIANGLEPHI
#define DEF_TRIANGLEPHI

#include "Triangle.hpp"

class TrianglePhi: public Triangle<std::complex<double> >{
	public:
		TrianglePhi(System const& s, double const& phi);
		~TrianglePhi() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		double phi_;

		void compute_H();
		void display_results();
		
		Matrix<double> set_ab();
		unsigned int match_pos_in_ab(Vector<double> const& x) const { (void)(x); return 0; }
};
#endif
