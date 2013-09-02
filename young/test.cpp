#include "YoungTableau.hpp"

#include<iostream>

int main(){
	unsigned int N(3);
	std::vector<unsigned int> a;
	a.push_back(1);
	YoungTableau yt_a(a,N);

	std::vector<unsigned int> b;
	b.push_back(1);
	YoungTableau yt_b(b,N);

	std::vector<unsigned int> c;
	c.push_back(1);
	YoungTableau yt_c(c,N);

	std::vector<YoungTableau> result(yt_a*yt_b*yt_c);
	std::cout<<yt_a<<" x "<<std::endl;
	std::cout<<yt_b<<" x "<<std::endl;
	std::cout<<yt_c<<" = "<<std::endl;
	std::cout<<yt_a*yt_b*yt_c<<std::endl;
}

