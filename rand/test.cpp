#include "Rand.hpp"

int main(){
	Rand rdm(1e4);
	
	std::cout<<rdm.get()<<std::endl;
	std::cout<<rdm.get(10)<<std::endl;
	std::cout<<rdm.get(2)<<std::endl;
}
