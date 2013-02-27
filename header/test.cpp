#include "Header.hpp"
#include "Read.hpp"
#include "Write.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"
#include <complex>

void write_bin(){
	std::cout<<"écriture d'un fichier binaire"<<std::endl;
	Header h("ce quatrième fichier contient");
	Matrice<double> a(6.3,6);
	h.add("m1",a);
	Write write("data-4.bin");
	h.write(write);
	write<<a;
}

void read_bin(){
	std::cout<<"lecture d'un fichier binaire"<<std::endl;
	Read read("data");
	Header h(read);
	h.show();
	double a(0);
	Matrice<double> m;
	read>>a>>m;
	std::cout<<a<<std::endl;
	std::cout<<m<<std::endl;
}

int main(){
	write_bin();
	//read_bin();
}

