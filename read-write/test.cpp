#include "Read.hpp"
#include "Write.hpp"
#include "Matrix.hpp"

void write_bin(){
	std::cout<<"écriture d'un fichier binaire"<<std::endl;
	Write write("data.bin");
	double a(12.3);
	std::string str("ci dessus, une matrice 5x3 valant 2.5 \n\nci dessous, une matrice complex 3x2 valant 3.2+0.2i\n");
	Matrix<double> M(5,3,2.5);
	std::complex<double> c(3.1,0.20);
	Matrix<std::complex<double> > C(2,3,c);
	Matrix<std::string> tab_str(5,2,"bla");
	tab_str(2,1) = "nouvelle string";	
	write<<M<<str<<C<<a<<c<<tab_str;
}

void read_bin(){
	std::cout<<"lecture d'un fichier binaire"<<std::endl;
	Read read_1("data.bin");
	Read read_2("data.bin");
	double a(0);
	std::string str1("");
	std::string str2("");
	Matrix<double> M1(5,4);
	Matrix<double> M2;
	std::complex<double> c(0.0);
	Matrix<std::complex<double> > C1(2,3);
	Matrix<std::complex<double> > C2;
	Matrix<std::string> tab_str1(5,2);
	Matrix<std::string> tab_str2;
	read_1>>M1>>str1>>C1>>a>>c>>tab_str1;

	std::cout<<M1<<std::endl;
	std::cout<<str1<<std::endl;
	std::cout<<C1<<std::endl;
	std::cout<<a<<" "<<c<<std::endl;
	std::cout<<tab_str1<<std::endl;

	std::cout<<"deuxième lecture sans mettre les dimensions aux matrix"<<std::endl;

	read_2>>M2>>str2>>C2>>a>>c>>tab_str2;
	std::cout<<M2<<std::endl;
	std::cout<<str2<<std::endl;
	std::cout<<C2<<std::endl;
	std::cout<<a<<" "<<c<<std::endl;
	std::cout<<tab_str2<<std::endl;
}

////void write_txt(){
	////Write write_single("single",false);
	////Write write_mat("matrice",false);
	////Write write_matc("matrice-complex",false);
	////Write write_arr2("array2d",false);
////
	////double a(12.3);
	////Matrice<double> M(5,2.5);
	////std::complex<double> c(3.1,0.20);
	////Matrice<std::complex<double> > C(2,c);
	////Array2D<int> A(5,2,-3);
////
	////write_single<<a<<c<<Write::endl;
	////write_mat <<M<<Write::endl;
	////write_matc<<C<<Write::endl;
	////write_arr2<<A<<Write::endl;
////
	////unsigned int po(18);
	////Write w;
	////w<<po;
	////w.open("post-ouverture",false);
	////w<<po;
////}

//void read_txt(){
	//Read read_single("single.dat");
	//Read read_mat("matrice.dat");
	//Read read_matc("matrice-complex.dat");
	//Read read_arr2("array2d.dat");
//
	//double a(0);
	//Matrice<double> M1(5);
	//std::complex<double> c(0.0);
	//Matrice<std::complex<double> > C1(2);
	//Array2D<int> A1(5,2);
//
	//read_single >> a >> c;
	//std::cout<<a<<" "<<c<<std::endl;
//
	//read_mat >> M1;
	//std::cout<<M1<<std::endl;
//
	//read_matc >> C1;
	//std::cout<<C1<<std::endl;
//
	//read_arr2 >> A1;
	//std::cout<<A1<<std::endl;
//
	//unsigned int po;
	//Read r;
	//r>>po;
	//std::cout<<"po="<<po<<std::endl;
	//r.open("post-ouverture",false);
	//r>>po;
	//std::cout<<"po="<<po<<std::endl;
//}

int main(){
	write_bin();
	read_bin();
	//write_txt();
	//read_txt();

}

