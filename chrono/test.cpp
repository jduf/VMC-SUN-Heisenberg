#include "Chrono.hpp"

int main(){
	Chrono t;
	t.tic();
	double r(0.0);
	for(int i(0);i<100;i++){
		std::cout<<t<<std::endl;
		for(int j(0);j<10000000;j++){
			r += i*j;
		}
		if(t.time_limit_reached(10)){ 
			std::cout<<"time limit reached"<<std::endl;
			i = 100000;
		}
	}

	t.tac();
	std::cout<<"final time "<<t<<std::endl;
}

