#ifndef DEF_HONEYCOMBSU4
#define DEF_HONEYCOMBSU4

#include "Honeycomb.hpp"

class HoneycombSU4: public Honeycomb<double>{
	public:
		HoneycombSU4(Container const& param);
		~HoneycombSU4();

		void create(double x);
		void study();
		void save();

	protected:
		void compute_T();
		void compute_P();
};
#endif

