#include "Read.hpp"
#include "Write.hpp"
#include "Matrix.hpp"
#include <complex>

void write_bin(){
	std::cout<<"écriture d'un fichier binaire"<<std::endl;
	double a(6.5);
	Matrix<double> A(10,2,2.2);
	std::complex<double> c(7.5,1.5);
	Matrix<std::complex<double> > C(3,8,c);

	Write write("data-3.jdbin");
	write("a",a);	
	write("A",A);	
	write("C",C);	
}

void read_bin(){
	std::cout<<"lecture d'un fichier binaire"<<std::endl;
	double a;
	Matrix<double> A;
	std::complex<double> c;
	Matrix<std::complex<double> > C;

	Read read("data-3.jdbin");
	std::cout<<read.get_header();
	read>>a>>A>>C;
	std::cout<<a<<std::endl;
	std::cout<<A<<std::endl;
	std::cout<<C<<std::endl;
}

int main(){
	write_bin();
	read_bin();
}

