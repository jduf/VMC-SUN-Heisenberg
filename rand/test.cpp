#include "Rand.hpp"
#include "Time.hpp"
#include<omp.h>

int main(){
	Rand rdm(1e4);
	
	std::cout<<rdm.get()<<std::endl;
	std::cout<<rdm.get(10)<<std::endl;
	std::cout<<rdm.get(2)<<std::endl;

	Time t;
	Rand rnd1(10);
	for(unsigned int i(0);i<1e5;i++){
		rnd1.get();
	}
	std::cout<<t.elapsed()<<std::endl;

	t.reset();
	Rand rnd2(1e6);
	for(unsigned int i(0);i<1e5;i++){
		rnd2.get();
	}
	std::cout<<t.elapsed()<<std::endl;

#pragma omp parallel for
	for(unsigned int i=0;i<4;i++){
		unsigned int thread(omp_get_thread_num());
		Rand rnd(1e4,thread);
		for(unsigned int j(0);j<5;j++){
			std::cout<<"thread="<<thread<<" rand="<<rnd.get(100)<<std::endl;;
		}
	}

	std::cout<<"new test"<<std::endl;
	Rand rnd3(1e4);
#pragma omp parallel for
	for(unsigned int i=0;i<20;i++){
		for(unsigned int j(0);j<10;j++){
			std::cout<<"thread="<<omp_get_thread_num()<<" rand="<<rnd3.get(100)<<std::endl;;
		}
	}
}
