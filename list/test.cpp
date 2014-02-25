#include "List.hpp"

int main(){
	List<double> a;
	for(unsigned int i(0);i<10;i++){
		a.append(i);
	}
	std::cout<<a<<std::endl;
	a.remove(4);
	std::cout<<a<<std::endl;
	a.remove(a.size()-1);
	std::cout<<a<<std::endl;
	a.remove(0);
	std::cout<<a<<std::endl;
	List<double> b;
	b.append(3);
	std::cout<<b<<std::endl;
	b.remove(0);
	std::cout<<b<<std::endl;
}
