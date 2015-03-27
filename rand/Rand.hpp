#ifndef DEF_RANDOM
#define DEF_RANDOM

#include <random>

/*{Description*/
/*!Random number generator that uses <random> of c++11. 
 * \warning if declared outside an parallel region but used by different
 * threads, to same random number may be used by different thread (openmp)
 */
/*}*/
template<typename Type>
class Rand{
	public:
		/*!Constructor that takes the range of generated random numbers*/
		Rand(Type const& min_inclusive, Type const& max_exclusive);
		/*!Destructor*/
		~Rand(){}
		/*!Gives the next random unsigned int strictly smaller than max*/
		Type get() { return dist_(mt_); }

	private:
		/*!Forbids default*/
		Rand();
		/*!Forbids copy*/
		Rand(Rand const&);

		mutable std::mt19937_64 mt_;
		mutable typename std::conditional< 
			std::is_integral<Type>::value,
			std::uniform_int_distribution<Type>,
			std::uniform_real_distribution<Type> >
				::type dist_;
};

template<typename Type>
Rand<Type>::Rand(Type const& min_inclusive, Type const& max_exclusive):
	mt_(std::random_device()()),
	dist_(min_inclusive,max_exclusive)
{}
#endif
