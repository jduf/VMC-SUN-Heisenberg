#ifndef DEF_PRIME
#define DEF_PRIME

#include<vector>
#include <iostream>

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

Prime::Prime(unsigned int N):
	p(new unsigned int[N]),
	N(N)
{
	p[0] = 2;
	unsigned int i(3),j(0);
	unsigned int new_prime(0);
	bool prime(true);
	while(new_prime<N-1){
		j = 0;
		while(p[j]*p[j] <= i && j<new_prime  && prime){
			if(i % p[j] == 0 ){prime=false;}
		j++;
		}

		if(prime){
			new_prime++;
			p[new_prime] = i;
		}
		prime=true;
		i++;
	}
}

Prime::~Prime(){
	delete[] p;
}

std::vector<unsigned int> Prime::pnd(unsigned int a) const{
	std::vector<unsigned int> out;
	out.push_back(1);
	if(p[N-1]*p[N-1]<a){
		std::cerr<<"Prime : not enough prime numbers to achive the prime number decomposition"<<std::endl;
	} else {
		unsigned int i(0);
		while(a>2){
			if(a % p[i] == 0){
				out.push_back(p[i]);
				a /= p[i];
			} else {
				i++;
			}
		}
	}
	return out;
}
#endif
