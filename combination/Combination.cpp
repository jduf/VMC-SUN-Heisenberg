#include "Combination.hpp"

Combination::Combination():
	k_(0),
	n_(0),
	comb_(NULL)
{}

void Combination::set(unsigned int const& k, unsigned int const& n, Vector<unsigned int>& comb){
	k_ = k;
	n_ = n;
	comb.set(k_);
	comb_=comb.ptr();
	for(unsigned int i(0);i<k_;i++){ comb(i) = i; }
}

bool Combination::next(){
	if(comb_){
		unsigned int i = k_ - 1;
		++comb_[i];
		/*the condition here may be problematic. I use the fact that i is an
		 * unsigned int and therefore whe i-1<0 is give a number that is very
		 * likely bigger than k. The original condition was with int i and i>=0*/
		while ((i-1 < k_) && (comb_[i] >= n_-k_+1+i)) {
			--i;
			++comb_[i];
		}
		/* Combination (n-k, n-k+1, ..., n) reached */
		/* No more combinations can be generated */
		if (comb_[0] > n_ - k_){ return false; }
		/* comb now looks like (..., x, n, n, n, ..., n).
		   Turn it into (..., x, x + 1, x + 2, ...) */
		for (i = i + 1; i < k_; ++i){ comb_[i] = comb_[i - 1] + 1;}

		return true;
	} else {
		std::cerr<<"bool Combination::next() : no vector to output the result"<<std::endl;
		return false;
	}
}

unsigned int Combination::number() const{
	unsigned int p(1);
	for(unsigned int i(0);i<n_-k_;i++){ p *= n_-i; }
	for(unsigned int i(1);i<n_-k_+1;i++){ p /= i; }
	return p;
}
