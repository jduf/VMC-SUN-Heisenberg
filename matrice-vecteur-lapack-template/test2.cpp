/**
@file test.cpp 
*/

#include "SquareMatrix.hpp"
#include "Vector.hpp"

#include <complex>

int main(){
	/*operateurs*/
	/*{*/
	unsigned int N(2);
	SquareMatrix<double> m1(N);
	SquareMatrix<double> m2(N,2);
	//
	m1(0,0) = 1;
	m1(0,1) = 2;
	m1(1,0) = 3;
	m1(1,1) = 4;
	m2(0,1) = 5;
	m2(1,0) = 6;
	m2(1,1) = 7;

	m2.test();
	//std::cout<<"m1"<<std::endl;
	//std::cout<<m1<<std::endl;
	//std::cout<<"m2"<<std::endl;
	//std::cout<<m2<<std::endl;
	//std::cout<<"m1*m2"<<std::endl;
	//std::cout<<m1*m2<<std::endl;
	//m1 *= m2;
	//std::cout<<m1<<std::endl;
	//
	//std::cout<<"m2-(m1*m2)"<<std::endl;
	//std::cout<<m2-m1<<std::endl;
	//m2-=m1;
	//std::cout<<"(m2-(m1*m2))-(m1*m2)"<<std::endl;
	//std::cout<<m2-m1<<std::endl;
	 /*}*/
}

