#ifndef RANMAR
#define RANMAR

#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <ctime>
#include <omp.h>

class Rand{
	public:

		// destructor
		~Rand();

		// these are the constructors

		// first constructor initialises new generator with length, seeds ij and kl
		Rand(int len, int ij, int kl);

		// next constructor initialises from save file fn
		//Rand(std::string fn, bool auto_save = true);

		// last constructor initialises generator with length len from seeds determined from system time
		Rand(int len);

		// constructor initialises generator with length len from seeds determined from system time but different on any given thread in the system
		Rand(int len, int threadnumber);

		// gives the next random number from the runvector
		inline float get() {
			pos++;
			if (pos == len){ ranvec();}
			return rvec[pos];
		}

		inline unsigned int get(unsigned int max) {
			pos++;
			if (pos == len){ ranvec();}
			return std::floor(rvec[pos] * max);
		}

		bool is_null();

		//void set_save(std::string fn) { Rand::fn = fn; auto_save = true; }

		// saves the state into a file
		//void save_state(std::string fn);

		// loades the state from a file, returns true if succeeded, false otherwise
		//bool load_state(std::string fn);

		// explicitly ask to shutdown
		void free(); 

	private:

		Rand() {};

		// ranvec calculates the next vector of random numbers
		void ranvec();

		void stop(std::string msg);

		// get seeds from system time
		void sow(int *, int*);

		float u[97], c, cd, cm;
		int i97, j97;

		int len; // the length of random vector
		int pos;
		float* rvec;

		//bool auto_save;
		//std::string fn;
};
#endif 
