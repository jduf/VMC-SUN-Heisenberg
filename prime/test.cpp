#include "Prime.hpp"

#include <iostream>
#include <vector>

int main(){
	unsigned int N(50);
	std::cout<<"N=";
	std::cin>>N;
	Prime P(10);
	std::vector<unsigned int> pnd(P.pnd(N));
	for(unsigned int i(0); i<pnd.size();i++){
		std::cout<<pnd[i]<<" ";
	}
	std::cout<<std::endl;
}



