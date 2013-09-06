#include "YoungTableau.hpp"

#include<iostream>

int main(){
	unsigned int N(0);
	unsigned int n_tab(0);
	do{
		std::cout<<"How many tableaus do you want to multiply ? : ";
		std::cin>>n_tab;
		std::cin.clear();
		std::cin.ignore(100,'\n');
	} while ( !std::cin );
	do{
		std::cout<<"SU(?) : ";
		std::cin>>N;
		std::cin.clear();
		std::cin.ignore(100,'\n');
	} while ( !std::cin || N<2);
	if (n_tab>1){
		std::vector<YoungTableau> tab_yt;
		for(unsigned int i(0);i<n_tab;i++){
			tab_yt.push_back(YoungTableau (N));
			std::cout<<tab_yt[i]<<std::endl;
		}

		std::vector<YoungTableau> result(tab_yt[0]*tab_yt[1]);
		for(unsigned int i(2);i<n_tab;i++){
			std::vector<YoungTableau> tmp(result*tab_yt[i]);
			for(unsigned int j(0);j<tmp.size();j++){
				result.push_back(tmp[j]);
			}
		}
		std::cout<<"*** will multiply the input tableaus ***"<<std::endl;
		for(unsigned int i(0);i<n_tab-1;i++){
			std::cout<<tab_yt[i]<<" x "<<std::endl;
		}
		std::cout<<tab_yt[n_tab-1]<<" = "<<std::endl;
		std::cout<<result<<std::endl;
	} else {
		std::cout<<YoungTableau(N)<<std::endl;
	}

}

