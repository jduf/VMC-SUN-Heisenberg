#ifndef DEF_PRIME
#define DEF_PRIME

#include <iostream>
#include <vector>

class Prime{
	public:
		Prime(unsigned int N);
		~Prime();

		inline unsigned int operator[](unsigned int const& i) const {return p[i];};
		std::vector<unsigned int> pnd(unsigned int a) const; //prime number decomposition

	private:
		Prime(Prime const& P);
		unsigned int *p;
		unsigned int N;
};

#endif
