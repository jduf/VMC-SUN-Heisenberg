#ifndef RANMAR
#define RANMAR

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <omp.h>

/*{Description*/
/*!
 * This random number generator originally appeared in "Toward a Universal
 * Random Number Generator" by George Marsaglia and Arif Zaman.
 * Florida State University Report: FSU-SCRI-87-50 (1987)
 * 
 * It was later modified by F. James and published in "A Review of Pseudo-
 * random Number Generators"
 * 
 * Converted from FORTRAN to C by Phil Linttell, James F. Hickling
 * Management Consultants Ltd, Aug. 14, 1989.
 * 
 * THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
 * (However, a newly discovered technique can yield
 * a period of 10^600. But that is still in the development stage.)
 * 
 * It passes ALL of the tests for random number generators and has a period
 * of 2^144, is completely portable (gives bit identical results on all
 * machines with at least 24-bit mantissas in the floating point
 * representation).
 * 
 * The algorithm is a combination of a Fibonacci sequence (with lags of 97
 * and 33, and operation "subtraction plus one, modulo one") and an
 * "arithmetic sequence" (using subtraction).
 * 
 * On a Vax 11/780, this random number generator can produce a number in
 * 13 microseconds.
 * 
 * The initialization routine for the random number generator RANMAR()
 * NOTE: The seed variables can have values between:   
 * 0 <= IJ <= 31328
 * 0 <= KL <= 30081
 * The random number sequences created by these two seeds are of sufficient
 * length to complete an entire calculation with. For example, if several
 * different groups are working on different parts of the same calculation,
 * each group could be assigned its own IJ seed. This would leave each group
 * with 30000 choices for the second seed. That is to say, this random
 * number generator can create 900 million different subsequences -- with
 * each subsequence having a length of approximately 10^30.
 * 
 * Use IJ = 1802 & KL = 9373 to test the random number generator. The
 * subroutine RANMAR should be used to generate 20000 random numbers.
 * Then display the next six random numbers generated multiplied by 4096*4096
 * If the random number generator is working properly, the random numbers
 * should be:
 * 6533892.0  14220222.0   7275067.0
 * 6172232.0   8354498.0  10633180.0
}*/
class Rand{
	public:
		/*!Constructor initialises new generator with length, seeds ij and kl*/
		Rand(int len, int ij, int kl);
		/*!Constructor initialises generator with length len from seed given by an other Rand*/
		Rand(int len, Rand& seed);
		/*!Constructor initialises generator with length len from seeds determined from system time*/
		Rand(int len);
		/*!Constructor initialises generator with length len from seeds determined from system time but different on any given thread in the systemr*/
		Rand(int len, int threadnumber);
		/*!Destructor*/
		~Rand();

		/*!Gives the next random float*/
		float get() {
			pos_++;
			if (pos_ == len_){ ranvec();}
			return rvec_[pos_];
		}

		/*!Gives the next random unsigned int strictly smaller than max*/
		unsigned int get(unsigned int max){
			pos_++;
			if (pos_ == len_){ ranvec();}
			return std::floor(rvec_[pos_] * max);
		}

	private:
		/*!Forbids default*/
		Rand(){};

		/*!Initializes the parameters*/
		void init(int seed1, int seed2);
		/*!Calculates the next vector of random numbers*/
		void ranvec();
		/*!Calls exit(0) and print msg*/
		void stop(std::string msg);
		/*{Description*/
		/*!The sow() procedure calculates two seeds for use with the random
		 * number generator from the system clock.  I decided how to do this
		 * myself, and I am sure that there must be better ways to select
		 * seeds; hopefully, however, this is good enough.  The first seed is
		 * calculated from the values for second, minute, hour, and year-day;
		 * weighted with the second most significant and year-day least
		 * significant.  The second seed weights the values in reverse.  */
		/*}*/
		void sow(int& seed1, int& seed2);

		float u_[97], c_, cd_, cm_;
		int i97_, j97_;

		int len_;	//!<the length of random vector
		float*rvec_;//!<array containing the random numbers
		int pos_;	//!<current position in the array
};
#endif 
