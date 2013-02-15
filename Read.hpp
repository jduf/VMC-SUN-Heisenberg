#ifndef DEF_READ
#define DEF_READ

#include"Matrice.hpp"
#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>


class Read{
	public:
		Read(std::string filename, bool binary=true);
		~Read();

		template<typename T>
			Read& operator>>(T& t);
		Read& operator>>(Matrice<double>& m);

	private:
		/*!forbids default constructor*/
		Read();
		/*!forbids copy constructor*/
		Read(Read const& R);
		/*!forbids assertion operator*/
		Read& operator=(Read const&);

		void open_binary(std::string filename);
		void read_binary_matrix(Matrice<double>& m);

		void open_txt(std::string filename);
		void read_txt_matrix(Matrice<double>& m);

		FILE *bfile;
		std::ifstream tfile;
		bool locked;
		bool binary;
};

template<typename T>
Read& Read::operator>>(T& t){
	if(binary){ fread(&t,sizeof(t),1,bfile); }
	else { tfile >> t; }
	return (*this);
}
#endif

