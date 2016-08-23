#ifndef DEF_RANDOM
#define DEF_RANDOM

#include <random>
#include <cassert>

/*{Description*/
/*!Random number generator that uses <random> of c++11.
 * \warning if declared outside an parallel region but used by different
 * threads, to same random number may be used by different thread (openmp)
 */
/*}*/
template<typename Type>
class Rand{
	public:
		/*!Constructor that sets [min,max) for double and [min,max] for unsigned int*/
		Rand(Type const& min, Type const& max);
		/*!Default destructor*/
		~Rand() = default;
		/*{Forbidden*/
		Rand() = delete;
		Rand(Rand const&) = delete;
		Rand(Rand&&) = delete;
		Rand& operator=(Rand) = delete;
		/*}*/

		/*!Gives the next random number*/
		Type get() const { return dist_(mt_); }

	private:
		mutable std::mt19937_64 mt_;
		mutable typename std::conditional<
			std::is_integral<Type>::value,
			std::uniform_int_distribution<Type>,
			std::uniform_real_distribution<Type> >
				::type dist_;
};

template<typename Type>
Rand<Type>::Rand(Type const& min, Type const& max):
	mt_(std::random_device()()),
	dist_(min,max)
{}

template<typename Type>
class RandArray{
	public:
		/*!Default constructor*/
		RandArray();
		/*!Constructor allocates memory for size_ rnd*/
		RandArray(unsigned int const& size);
		/*!Dedfault destructor*/
		~RandArray();
		/*{Forbidden*/
		RandArray(RandArray const&) = delete;
		RandArray(RandArray&&) = delete;
		RandArray& operator=(RandArray) = delete;
		/*}*/

		void set(unsigned int const& size);
		/*!Set the ith rnd to generate numbers within [min,max) for double and
		 * [min,max] for unsigned int*/
		void set(unsigned int const& i, Type const& min, Type const& max);
		/*!Gives the next random number for the ith index*/
		Type get(unsigned int const& i) const { assert(i<size_); return rnd_[i]->get(); }

	private:
		unsigned int size_;
		Rand<Type>** rnd_;
};

template<typename Type>
RandArray<Type>::RandArray():
	size_(0),
	rnd_(NULL)
{}

template<typename Type>
RandArray<Type>::RandArray(unsigned int const& size):
	size_(size),
	rnd_(new Rand<Type>*[size])
{
	for(unsigned int i(0);i<size_;i++){ rnd_[i] = NULL; }
}

template<typename Type>
RandArray<Type>::~RandArray(){
	for(unsigned int i(0);i<size_;i++){ if(rnd_[i]){ delete rnd_[i]; } }
	delete[] rnd_;
	rnd_ = NULL;
}

template<typename Type>
void RandArray<Type>::set(unsigned int const& size){
	if(rnd_){
		for(unsigned int i(0);i<size_;i++){
			delete rnd_[i];
			rnd_[i] = NULL;
		}
	}
	size_ = size;
	rnd_ = new Rand<Type>*[size];
	for(unsigned int i(0);i<size_;i++){ rnd_[i] = NULL; }
}

template<typename Type>
void RandArray<Type>::set(unsigned int const& i, Type const& min, Type const& max){
	assert(i<size_);
	if(rnd_[i]){ delete rnd_[i]; }
	rnd_[i] = new Rand<Type>(min,max);
}
#endif
