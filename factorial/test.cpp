#include "Time.hpp"
#include <iostream>

int main(){
	Time t;
	unsigned int n(25);
	int f(1);
	for(unsigned int i(1);i<n;i++){
		f *= i;
		std::cout<<i<<" "<<f<<std::endl;
	}
	std::cout<<"final time "<<t.elapsed()<<std::endl;
}

