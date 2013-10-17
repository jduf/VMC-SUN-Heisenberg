#include "YoungTableau.hpp"

#include<iostream>


int main(){
	unsigned int N(3);
	unsigned int N_s(4);
	std::vector<unsigned int> a;
	a.push_back(2);
	a.push_back(1);
	YoungTableau yt_a(a,N);
	std::vector<YoungTableau> list;
	list.push_back(yt_a);
	for(unsigned int k(1); k<N_s;k++){
		std::vector<YoungTableau> tmp;
		for(unsigned int i(0);i<list.size();i++){
			std::vector<YoungTableau> new_list(list[i]*yt_a);
			for(unsigned int j(0); j<new_list.size();j++){
				tmp.push_back(new_list[j]);
			}
		}
		list = tmp;
		std::cout<<k<<std::endl;
	}
	unsigned int N_singulet(0);
	for(unsigned int i(0);i<list.size();i++){
		if(list[i].dimension()==1){ N_singulet++; }
	}
	std::cout<<N_singulet<<std::endl;
	//std::cout<<list<<std::endl;
}
