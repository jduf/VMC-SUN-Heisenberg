#ifndef DEF_HONEYCOMBSU4
#define DEF_HONEYCOMBSU4

#include "Honeycomb.hpp"

class HoneycombSU4: public Honeycomb<double>{
	public:
		HoneycombSU4(unsigned int N, unsigned int n, unsigned int m);
		~HoneycombSU4();

		void create(double x);
		void check();
		void study();

	protected:
		void compute_T();
};
#endif

