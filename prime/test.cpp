#include "Prime.hpp"

#include<cstdlib>
int main(){
	unsigned int N(50);
	Prime P(10);
	std::vector<unsigned int> pnd(P.pnd(N));
	for(unsigned int i(0); i<pnd.size();i++){
		std::cout<<pnd[i]<<" ";
	}
	std::cout<<std::endl;
}



