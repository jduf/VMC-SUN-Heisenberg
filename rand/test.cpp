#include "Rand.hpp"
#include "Chrono.hpp"

int main(){
	Rand rdm(1e4);
	
	std::cout<<rdm.get()<<std::endl;
	std::cout<<rdm.get(10)<<std::endl;
	std::cout<<rdm.get(2)<<std::endl;

	Chrono t1,t2;
	t1.tic();
	for(unsigned int i(0);i<1e5;i++){
		Rand rnd1(10);
		rnd1.get();
	}
	t1.tac();
	t2.tic();
	Rand rnd2(1e6);
	for(unsigned int i(0);i<1e5;i++){
		rnd2.get();
	}
	t2.tac();
	std::cout<<t1<<" "<<t2<<std::endl;
}
