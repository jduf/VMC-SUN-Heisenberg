#ifndef DEF_COMBINATION
#define DEF_COMBINATION

#include "Vector.hpp"

/*!taken from : 
 * http://compprog.wordpress.com/2007/10/17/generating-combinations-1/
 * 
 * next_comb(int comb[], int k, int n) Generates the next combination of n
 * elements as k after comb
 * 
 * comb => the previous combination ( use (0, 1, 2, ..., k) for first)
 * k => the size of the subsets to generate
 * n => the size of the original set
 * 
 * Returns: 1 if a valid combination was found
 * 0, otherwise
*/
class Combination{
	public: 
		Combination();
		void set(unsigned int const& k, unsigned int const& n,Vector<unsigned int>& comb);

		bool next();
		unsigned int number() const;

	private:
		unsigned int k_;
		unsigned int n_;
		unsigned int* comb_;
};
#endif
