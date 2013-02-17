#ifndef DEF_READ
#define DEF_READ

#include"Matrice.hpp"
#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include<complex>


class Read{
	public:
		Read(std::string filename, bool binary=true);
		~Read();

		template<typename T>
			Read& operator>>(T& t);
		template<typename M>
			Read& operator>>(Matrice<M>& m);

	private:
		/*!forbids default constructor*/
		Read();
		/*!forbids copy constructor*/
		Read(Read const& R);
		/*!forbids assertion operator*/
		Read& operator=(Read const&);

		void open_binary(std::string filename);
		void open_txt(std::string filename);

		template<typename M>
			void read_binary_matrix(Matrice<M>& m);
		template<typename M>
			void read_txt_matrix(Matrice<M>& m);

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

template<typename M>
Read& Read::operator>>(Matrice<M>& m){
	if(binary) { read_binary_matrix(m); }
	else { read_txt_matrix(m); }
	return (*this);
}

template<typename M>
void Read::read_binary_matrix(Matrice<M>& m){
	unsigned int N(0);
	fread(&N,sizeof(N),1,bfile);
	M tmp[N*N];
	fread(&tmp,sizeof(tmp),1,bfile);
	if(N != m.size()) {
		Matrice<M> mat_tmp(N);
		for(unsigned int i(0);i<N*N;i++){
			(mat_tmp.ptr())[i]=tmp[i];
		}
		m = mat_tmp;
	} else {
		for(unsigned int i(0);i<N*N;i++){
			(m.ptr())[i]=tmp[i];
		}
	}
}

template<typename M>
void Read::read_txt_matrix(Matrice<M>& m){
	if(m.size()!=0) {
		for(unsigned int i(0); i<m.size();i++){
			for(unsigned int j(0); j<m.size();j++){
				tfile >> m(i,j);
			}
		}
	} else {
		std::cerr<<"Read : to read a Matrice<M> you need to set the size of the input matric"<<std::endl;
	}
}
#endif

