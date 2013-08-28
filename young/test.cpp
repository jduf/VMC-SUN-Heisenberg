#include "YoungTableau.hpp"

#include<vector>


int main(){
	std::vector<unsigned int> a;
	a.push_back(2);
	a.push_back(1);
	YoungTableau yt_a(a);

	std::vector<unsigned int> b;
	b.push_back(2);
	b.push_back(1);
	YoungTableau yt_b(b);

	std::cout<<yt_a.dimension()<<std::endl;
	yt_a.print();
	std::cout<<yt_b.dimension()<<std::endl;
	yt_b.print();

	std::vector<YoungTableau> result(yt_a*yt_b);
	for(unsigned int i(0); i<result.size(); i++){
		result[i].print();
	}

}

