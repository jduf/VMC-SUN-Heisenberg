#ifndef DEF_WRITE
#define DEF_WRITE

#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include "Matrice.hpp"


class Write{
	public:
		Write(std::string filename, bool binary=true);
		~Write();


		template<typename T>
			Write& operator<<(T const& t);	
		Write& operator<<(Matrice<double> const& m);
		
	private:
		/*!forbids default constructor*/
		Write();
		/*!forbids copy constructor*/
		Write(Write const& s);
		/*!forbids assertion operator*/
		Write& operator=(Write const&);

		void open_binary(std::string filename);
		void write_binary_matrix(Matrice<double> const& m);

		void open_txt(std::string filename);
		void write_txt_matrix(Matrice<double> const& m);

		FILE *bfile;
		std::ofstream tfile;
		bool binary;
		bool locked;

};


template<typename T>
Write& Write::operator<<(T const& t){
	if(binary){
		fwrite(&t, sizeof(t), 1 ,bfile);
		fflush(bfile);
	} else { tfile<< t; }
	return (*this);
}

#endif

