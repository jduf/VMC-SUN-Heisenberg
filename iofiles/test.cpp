#include "IOFiles.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

void write_bin(){
	std::cout<<"écriture d'un fichier binaire"<<std::endl;
	Matrix<double> M(5,3,2.5);
	std::complex<double> c(1.2,2);
	Vector<std::complex<double> > v(3,c);
	double a(12.2);
	std::string slt("salut",false);
	IOFiles write("data.jdbin",true,false);
	write<<a<<M<<c<<slt;
	write.write("Vector",v);
}

void read_bin(){
	std::cout<<"lecture d'un fichier binaire"<<std::endl;
	IOFiles read("data.jdbin",false,false);
	Matrix<double> M;
	Vector<std::complex<double> > v;
	std::complex<double> c;
	double a;
	std::string slt;
	read>>a>>M>>c>>slt>>v;
	//std::cout<<read.get_header()<<std::endl;
	//std::cout<<a<<std::endl;
	//std::cout<<M<<std::endl;
	//std::cout<<c<<std::endl;
	//std::cout<<slt<<std::endl;
	//std::cout<<v<<std::endl;
}

void write_txt(){
	std::cout<<"écriture d'un fichier text"<<std::endl;
	IOFiles write_single("single.dat",true,false);
	IOFiles write_mat("matrix.dat",true,false);
	IOFiles write_vec_complex("vector-complex.dat",true,false);

	double a(12.3);
	Matrix<double> M(5,5,2.5);
	std::complex<double> c(3.1,0.20);
	Vector<std::complex<double> > C(6,c);

	write_single<<a<<" "<<c<<IOFiles::endl;
	write_mat <<M<<IOFiles::endl;
	write_vec_complex<<C<<IOFiles::endl;
}

void read_txt(){
	std::cout<<"lecture d'un fichier text"<<std::endl;
	IOFiles read_single("single.dat",false,false);
	IOFiles read_mat("matrix.dat",false,false);
	IOFiles read_vec_complex("vector-complex.dat",false,false);

	double a(0);
	Matrix<double> M1(5,5);
	std::complex<double> c(0.0);
	Vector<std::complex<double> > C1(6);

	read_single >> a >> c;
	std::cout<<a<<" "<<c<<std::endl;

	read_mat >> M1;
	std::cout<<M1<<std::endl;

	read_vec_complex >> C1;
	std::cout<<C1<<std::endl;
}

void append_txt(){
	std::cout<<"addition à un fichier text"<<std::endl;
	IOFiles write_single("single.dat",true,true);
	IOFiles write_mat("matrix.dat",true,true);
	IOFiles write_vec_complex("vector-complex.dat",true,true);

	write_single<<"extra"<<IOFiles::endl;
	write_mat<<"extra"<<IOFiles::endl;
	write_vec_complex<<"extra"<<IOFiles::endl;
}

int main(){
	//write_bin();
	//read_bin();
	write_txt();
	read_txt();
	append_txt();
}
