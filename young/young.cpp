#include "YoungTableau.hpp"

#include<vector>

int main(){
	unsigned int N(3);
	YoungTableau yt_a(N);
	YoungTableau yt_b(N);
	//std::vector<unsigned int> a;
	//a.push_back(2);
	//a.push_back(1);
	//YoungTableau yt_a(a,N);
//
	//std::vector<unsigned int> b;
	//b.push_back(2);
	//b.push_back(1);
	//YoungTableau yt_b(b,N);
//
	std::vector<YoungTableau> result(yt_a*yt_b);
	yt_a.print();
	std::cout<<" x "<<std::endl;
	yt_b.print();
	std::cout<<" = "<<std::endl;
	for(unsigned int i(0); i<result.size(); i++){
		result[i].print();
		std::cout<<" + "<<std::endl;
	}
}

