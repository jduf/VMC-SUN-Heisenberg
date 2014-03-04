#ifndef DEF_HONEYCOMBSU3
#define DEF_HONEYCOMBSU3

#include "Honeycomb.hpp"

class HoneycombSU3: public Honeycomb<double>{
	public:
		HoneycombSU3(unsigned int N, unsigned int n, unsigned int m);
		~HoneycombSU3();

		void create(double x);
		void check();
		void study();

	protected:
		void compute_T();
};
#endif

