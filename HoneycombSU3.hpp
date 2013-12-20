#ifndef DEF_HONEYCOMBSU3
#define DEF_HONEYCOMBSU3

#include "Honeycomb.hpp"

class HoneycombSU3: public Honeycomb<double>{
	public:
		HoneycombSU3(Parseur& P);
		~HoneycombSU3();

	protected:
		void compute_T();
		void save();
};
#endif

