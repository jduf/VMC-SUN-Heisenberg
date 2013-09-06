#include "YoungTableau.hpp"

#include<iostream>

int main(){
	unsigned int N(3);
	std::vector<unsigned int> a;
	a.push_back(4);
	a.push_back(2);
	YoungTableau yt_a(a,N);
	yt_a.test();
}
