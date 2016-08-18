#include "Vector.hpp"
#include "Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	unsigned int n(P.get<unsigned int>("n"));
	Vector<unsigned int> tmp;

	tmp = P.get<std::vector<unsigned int> >("p0");
	Matrix<unsigned int> p0(2,n);
	for(unsigned int i(0);i<n;i++){
		p0(0,i) = i+1; 
		p0(1,i) = i+1; 
	}
	for(unsigned int i(0);i<tmp.size()-1;i++){ p0(1,tmp(i)-1) = tmp(i+1); }
	p0(1,tmp.back()-1) = tmp(0);

	tmp = P.get<std::vector<unsigned int> >("p1");
	Matrix<unsigned int> p1(2,n);
	for(unsigned int i(0);i<n;i++){
		p1(0,i) = i+1; 
		p1(1,i) = i+1; 
	}
	for(unsigned int i(0);i<tmp.size()-1;i++){ p1(1,tmp(i)-1) = tmp(i+1); }
	p1(1,tmp.back()-1) = tmp(0);


	Matrix<unsigned int> p2(2,n);
	for(unsigned int i(0);i<n;i++){
		p2(0,i) = i+1; 
		p2(1,i) = p1(1,p0(1,i)-1); 
	}

	std::cout<<p0<<std::endl;
	std::cout<<std::endl;
	std::cout<<p1<<std::endl;
	std::cout<<std::endl;
	std::cout<<p2<<std::endl;
}
