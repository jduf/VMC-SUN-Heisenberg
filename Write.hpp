#ifndef DEF_WRITE
#define DEF_WRITE

#include "Matrice.hpp"
#include "Array2D.hpp"

#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include<complex>


class Write{
	public:
		Write(std::string filename, bool binary=true);
		Write();
		~Write();

		static std::string endl;

		template<typename T>
			Write& operator<<(T const& t);	
		template<typename M>
			Write& operator<<(Matrice<M> const& mat);
		template<typename A>
			Write& operator<<(Array2D<A> const& arr);
	
		void open(std::string filename, bool binary=true);

	private:
		/*!forbids copy constructor*/
		Write(Write const& s);
		/*!forbids assertion operator*/
		Write& operator=(Write const&);

		void open_binary();
		void open_txt();

		template<typename M>
			void write_binary_matrix(Matrice<M> const& mat);
		template<typename A>
			void write_binary_array2d(Array2D<A> const& arr);

		std::string filename;
		FILE *bfile;
		std::ofstream tfile;
		bool unlocked;
		bool binary;
};

template<typename T>
Write& Write::operator<<(T const& t){
	if(unlocked){
		if(binary){
			fwrite(&t, sizeof(t), 1 ,bfile);
			fflush(bfile);
		} else { tfile<< t; }
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename M>
Write& Write::operator<<(Matrice<M> const& mat){
	if(unlocked){
		if(binary) { write_binary_matrix(mat); }
		else { tfile<<mat; }
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename A>
Write& Write::operator<<(Array2D<A> const& arr){
	if(unlocked){
		if(binary) { write_binary_array2d(arr);}
		else { tfile<<arr; }
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename M>
void Write::write_binary_matrix(Matrice<M> const& mat){
	unsigned int N(mat.size());
	fwrite(&N, sizeof(N), 1 ,bfile);
	M tmp[N*N];
	for(unsigned int i(0);i<N*N;i++){
		tmp[i] = (mat.ptr())[i];
	}
	fwrite(&tmp, sizeof(tmp), 1 ,bfile);
	fflush(bfile);
}

template<typename A>
void Write::write_binary_array2d(Array2D<A> const& arr){
	unsigned int N_row(arr.row());
	unsigned int N_col(arr.col());

	fwrite(&N_row, sizeof(N_row), 1 ,bfile);
	fwrite(&N_col, sizeof(N_col), 1 ,bfile);
	A tmp[N_row*N_col];
	for(unsigned int i(0);i<N_row*N_col;i++){
		tmp[i] = (arr.ptr())[i];
	}
	fwrite(&tmp, sizeof(tmp), 1 ,bfile);
	fflush(bfile);
}
#endif

