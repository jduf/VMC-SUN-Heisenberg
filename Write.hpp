#ifndef DEF_WRITE
#define DEF_WRITE

#include "Matrice.hpp"
#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include<complex>


class Write{
	public:
		Write(std::string filename, bool binary=true);
		~Write();

		static std::string endl;

		template<typename T>
			Write& operator<<(T const& t);	
		template<typename M>
			Write& operator<<(Matrice<M> const& m);
		
	private:
		/*!forbids default constructor*/
		Write();
		/*!forbids copy constructor*/
		Write(Write const& s);
		/*!forbids assertion operator*/
		Write& operator=(Write const&);

		void open_binary(std::string filename);
		void open_txt(std::string filename);

		template<typename M>
			void write_binary_matrix(Matrice<M> const& m);
		template<typename M>
			void write_txt_matrix(Matrice<M> const& m);


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

template<typename M>
Write& Write::operator<<(Matrice<M> const& m){
	if(binary) { write_binary_matrix(m); }
	else { write_txt_matrix(m); }

	return (*this);
}

template<typename M>
void Write::write_binary_matrix(Matrice<M> const& m){
	unsigned int N(m.size());
	fwrite(&N, sizeof(N), 1 ,bfile);
	M tmp[N*N];
	for(unsigned int i(0);i<N*N;i++){
		tmp[i] = (m.ptr())[i];
	}
	fwrite(&tmp, sizeof(tmp), 1 ,bfile);
	fflush(bfile);
}

template<typename M>
void Write::write_txt_matrix(Matrice<M> const& m){
	for(unsigned int i(0); i<m.size();i++){
		for(unsigned int j(0); j<m.size();j++){
			tfile << m(i,j)<<" ";
		}
		tfile<<std::endl;
	}
}
#endif

