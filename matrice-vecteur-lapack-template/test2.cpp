/**
@file test.cpp 
*/

#include "Matrix.hpp"
#include "SquareMatrix.hpp"

#include <complex>

int main(){
	/*operateurs*/
	/*{*/
	unsigned int N1_row(3);
	unsigned int N1_col(2);
	unsigned int N2_row(2);
	unsigned int N2_col(6);
	Matrix<double> m1(N1_row,N1_col);
	Matrix<double> m2(N2_row,N2_col,2);
	//
	m1(0,0) = 1;
	m1(0,1) = 2;
	m1(1,0) = 3;
	m1(1,1) = 4;
	m1(2,0) = 5;
	m1(2,1) = 6;

	std::cout<<"m1"<<std::endl;
	std::cout<<m1<<std::endl;
	std::cout<<"m2"<<std::endl;
	std::cout<<m2<<std::endl;
	std::cout<<"m1*m2"<<std::endl;
	m1 = m1*m2;
	std::cout<<m1<<std::endl;
	m2=m2.transpose();
	std::cout<<"m2.tanspose"<<std::endl;
	std::cout<<m2<<std::endl;
	std::cout<<"(m1*m2)*m2"<<std::endl;
	std::cout<<m1*m2<<std::endl;

	
	unsigned int N(3);
	SquareMatrix<double> m3(N,3.3);

	std::cout<<m3<<std::endl;
	m3 *= m3;
	std::cout<<"m3^2"<<std::endl;
	std::cout<<m3<<std::endl;

	//
	//std::cout<<"m2-(m1*m2)"<<std::endl;
	//std::cout<<m2-m1<<std::endl;
	//m2-=m1;
	//std::cout<<"(m2-(m1*m2))-(m1*m2)"<<std::endl;
	//std::cout<<m2-m1<<std::endl;
	 /*}*/
}

