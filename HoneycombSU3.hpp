#ifndef DEF_HONEYCOMBSU3
#define DEF_HONEYCOMBSU3

#include "Honeycomb.hpp"

class HoneycombSU3: public Honeycomb<double>{
	public:
		HoneycombSU3(Container const& param);
		~HoneycombSU3();

		void save();

	protected:
		void compute_T();
};
#endif

