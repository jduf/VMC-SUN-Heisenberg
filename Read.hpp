#ifndef DEF_READ
#define DEF_READ

#include"Matrice.hpp"
#include"Array2D.hpp"

#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include<complex>


class Read{
	public:
		Read(std::string filename, bool binary=true);
		Read();
		~Read();

		template<typename T>
			Read& operator>>(T& t);
		template<typename M>
			Read& operator>>(Matrice<M>& mat);
		template<typename A>
			Read& operator>>(Array2D<A>& arr);
		
		void open(std::string filename, bool binary=true);

	private:
		/*!forbids copy constructor*/
		Read(Read const& R);
		/*!forbids assertion operator*/
		Read& operator=(Read const&);

		void open_binary();
		void open_txt();

		template<typename M>
			void read_binary_matrix(Matrice<M>& mat);
		template<typename A>
			void read_binary_array2d(Array2D<A>& arr);

		std::string filename;
		FILE *bfile;
		std::ifstream tfile;
		bool unlocked;
		bool binary;
};

template<typename T>
Read& Read::operator>>(T& t){
	if(unlocked){
		if(binary){ fread(&t,sizeof(t),1,bfile); }
		else { tfile >> t; }
	} else {
		std::cerr<<"Read : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename M>
Read& Read::operator>>(Matrice<M>& mat){
	if(unlocked){
		if(binary) { read_binary_matrix(mat); }
		else { tfile>>mat; }
	} else {
		std::cerr<<"Read : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename A>
Read& Read::operator>>(Array2D<A>& arr){
	if(unlocked){
		if(binary) { read_binary_array2d(arr);}
		else { tfile>>arr; }
	} else {
		std::cerr<<"Read : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename M>
void Read::read_binary_matrix(Matrice<M>& mat){
	unsigned int N(0);
	fread(&N,sizeof(N),1,bfile);
	M tmp[N*N];
	fread(&tmp,sizeof(tmp),1,bfile);
	if(N != mat.size()) {
		Matrice<M> mat_tmp(N);
		for(unsigned int i(0);i<N*N;i++){
			(mat_tmp.ptr())[i]=tmp[i];
		}
		mat = mat_tmp;
	} else {
		for(unsigned int i(0);i<N*N;i++){
			(mat.ptr())[i]=tmp[i];
		}
	}
}

template<typename A>
void Read::read_binary_array2d(Array2D<A>& arr){
	unsigned int N_row(0);
	unsigned int N_col(0);
	fread(&N_row,sizeof(N_row),1,bfile);
	fread(&N_col,sizeof(N_col),1,bfile);
	A tmp[N_row*N_col];
	fread(&tmp,sizeof(tmp),1,bfile);
	if(N_row != arr.row() || N_col != arr.col()) {
		Array2D<A> arr_tmp(N_row,N_col);
		for(unsigned int i(0);i<N_row*N_col;i++){
			(arr_tmp.ptr())[i]=tmp[i];
		}
		arr = arr_tmp;
	} else {
		for(unsigned int i(0);i<N_row*N_col;i++){
			(arr.ptr())[i]=tmp[i];
		}
	}
}
#endif

