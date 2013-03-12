#include "Read.hpp"
#include "Write.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"
#include <complex>

void write_bin(){
	std::cout<<"écriture d'un fichier binaire"<<std::endl;
	double a(6.5);
	Matrice<double> A(10,2.2);
	std::complex<double> c(7.5,1.5);
	Matrice<std::complex<double> > C(8,c);

	Write write("data-3.jdbin");
	write("a",a);	
	//write("c",c);	
	write("A",A);	
	write("C",C);	
	write.linking();
}

void read_bin(){
	std::cout<<"lecture d'un fichier binaire"<<std::endl;
	double a;
	Matrice<double> A;
	std::complex<double> c;
	Matrice<std::complex<double> > C;

	Read read("data-3.jdbin");
	std::cout<<read.header();
	//read>>a>>c>>A>>C;
	read>>a>>A>>C;
	//read>>a>>c;
	std::cout<<a<<std::endl;
	std::cout<<A<<std::endl;
	//std::cout<<c<<std::endl;
	std::cout<<C<<std::endl;
}


int main(){
	write_bin();
	read_bin();
}

