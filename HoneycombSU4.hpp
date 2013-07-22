#ifndef DEF_HONEYCOMBSU4
#define DEF_HONEYCOMBSU4

#include "Honeycomb.hpp"

class HoneycombSU4: public Honeycomb<double>{
	public:
		HoneycombSU4(Parseur& P);
		~HoneycombSU4();

	protected:
		void compute_EVec();
		void save(std::string filename);
};
#endif

