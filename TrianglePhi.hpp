#ifndef DEF_TRIANGLEPHI
#define DEF_TRIANGLEPHI

#include "Triangle.hpp"

class TrianglePhi: public Triangle<std::complex<double> >{
	public:
		TrianglePhi(Parseur& P);
		~TrianglePhi();

	protected:
		double phi_;
		
		void compute_T();
		void compute_P();
		void compute_band_structure();
		void save();
};
#endif

