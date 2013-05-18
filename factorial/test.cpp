#include "Chrono.hpp"

int main(){
	Chrono t;
	t.tic();
	unsigned int n(20);
	unsigned long int f(1);
	for(int i(1);i<20;i++){
		f *= i;
		std::cout<<i<<" "<<f<<std::endl;
	}
	t.tac();
	std::cout<<"final time "<<t<<std::endl;

}

