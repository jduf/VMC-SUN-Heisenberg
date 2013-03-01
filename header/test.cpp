#include "Read.hpp"
#include "Write.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"
#include <complex>

void write_bin(){
	std::cout<<"écriture d'un fichier binaire"<<std::endl;
	double a(7.35);
	Matrice<double> A(6,6.2);
	std::complex<double> c(2.3,3.5);
	Matrice<std::complex<double> > C(2,c);

	Write write("data-1.jdbin",true);
	write("a",a);	
	write("c",c);	
	write("A",A);	
	write("C",C);	
	//write<<a<<C<<c<<A;
}

void read_bin(){
	std::cout<<"lecture d'un fichier binaire"<<std::endl;
	double a;
	Matrice<double> A;
	std::complex<double> c;
	Matrice<std::complex<double> > C;

	Read read("data-1.jdbin",true);
	std::cout<<read.header();
	read>>a>>c>>A>>C;
	std::cout<<a<<std::endl;
	std::cout<<A<<std::endl;
	std::cout<<c<<std::endl;
	std::cout<<C<<std::endl;
}


int main(){
	write_bin();
	read_bin();
}

