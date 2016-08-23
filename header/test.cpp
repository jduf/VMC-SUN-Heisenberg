#include "IOFiles.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

void write_bin(){
	std::cout<<"écriture d'un fichier binaire"<<std::endl;
	double a(6.5);
	Matrix<double> A(10,2,2.2);
	std::complex<double> c(7.5,1.5);
	Vector<std::complex<double> > C(3,c);

	IOFiles write("data.jdbin",true,false);
	write.write("a",a);	
	write.write("A",A);	
	write.write("C",C);	
}

void read_bin(){
	std::cout<<"lecture d'un fichier binaire"<<std::endl;
	double a;
	Matrix<double> A;
	std::complex<double> c;
	Vector<std::complex<double> > C;

	IOFiles read("data.jdbin",false,false);
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

