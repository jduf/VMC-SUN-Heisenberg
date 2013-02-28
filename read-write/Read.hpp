#ifndef DEF_READ
#define DEF_READ

#include "Matrice.hpp"
#include "Array2D.hpp"
#include "Header.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <complex>

/*!Class that allows to read datas easily from a file.
 *To be used with Write.hpp
 *Can be used to read a binary or a text file.
 *
 *When used with a text file
 * - the datas must be separated with a space
 * - only one matrix or array is safely read per file
 * - the size of the matrix or array has to be specified
 *
 *When used with a binary file
 * - many different kinds of datas can be extracted
*/
class Read{
	public:
		/*!Default constructor that needs a call of Read::open(std::string filename, bool binary)*/
		Read();
		/*!Opens a file named "filename", reads from the filename the type of file*/
		Read(std::string filename, bool header=false);
		/*!Closes the file*/
		~Read();

		/*!Stream operator that reads datas without formatting*/
		template<typename T>
			Read& operator>>(T& t);
		/*!Stream operator that reads matrices, uses Matrice<M>::operator>>*/
		template<typename M>
			Read& operator>>(Matrice<M>& mat);
		/*!Stream operator that reads arrays, uses Array2D<A>::operator>>*/
		template<typename A>
			Read& operator>>(Array2D<A>& arr);
		/*!Stream operator that reads strings*/
		Read& operator>>(std::string& s);

		/*!To be used with a default constructor : opens a file named "filename", reads from the filename the type of file*/
		void open(std::string filename, bool header=false);

		inline std::string header() const { return (h->get())->get(); };

	private:
		/*!Forbids copy constructor*/
		Read(Read const& R);
		/*!Forbids assertion operator*/
		Read& operator=(Read const&);

		/*!Subroutine needed to open a binary file*/
		void open_binary();
		/*!Subroutine needed to open a text file*/
		void open_txt();
		/*!Extract the header from the file and save it in h*/
		void read_header();

		/*!Subroutine that check if the correct extension is given*/
		bool is_binary(std::string f);

		/*!Subroutine needed to read a matrix from a binary file*/
		template<typename M>
			void read_binary_matrix(Matrice<M>& mat);
		/*!Subroutine needed to read a matrix from a binary file*/
		template<typename A>
			void read_binary_array2d(Array2D<A>& arr);

		std::string filename; //!< name of the file to read from
		FILE *bfile; //!< pointer on the binery file to read from
		std::ifstream tfile; //!< stream of the text file to read form
		Header *h; //!< pointer to a header (actually it will be a footer)
		bool unlocked; //!< true if the file is ready to be read from
		bool binary; //!< true if the file is binary
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

