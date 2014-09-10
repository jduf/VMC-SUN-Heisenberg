#ifndef RANDOM
#define RANDOM

#include <random>

class RandDouble{
	public:
		/*!Constructor initialises generator with length, seeds ij and kl*/
		RandDouble(double const& min_inclusive, double const& max_exclusive);
		/*!Destructor*/
		~RandDouble(){};
		/*!Gives the next random unsigned int strictly smaller than max*/
		double get(){ return dist_(mt_); }

	private:
		/*!Forbids default*/
		RandDouble();
		/*!Forbids copy*/
		RandDouble(RandDouble const&);

		std::mt19937_64 mt_;
		std::uniform_real_distribution<double> dist_;
};

class RandUnsignedInt{
	public:
		/*!Constructor initialises generator with length, seeds ij and kl*/
		RandUnsignedInt(unsigned int const& min_inclusive, unsigned int const& max_inclusive);
		/*!Destructor*/
		~RandUnsignedInt(){};
		/*!Gives the next random unsigned int strictly smaller than max*/
		unsigned int get(){ return dist_(mt_); }

	private:
		/*!Forbids default*/
		RandUnsignedInt();
		/*!Forbids copy*/
		RandUnsignedInt(RandUnsignedInt const&);

		std::mt19937_64 mt_;
		std::uniform_int_distribution<unsigned int> dist_;
};
#endif 
