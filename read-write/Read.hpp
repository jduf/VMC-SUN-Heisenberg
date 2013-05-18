#ifndef DEF_READ
#define DEF_READ

#include "Matrice.hpp"
#include "Array2D.hpp"
#include "Header.hpp"

#include <iostream>
#include <fstream>
#include <string>
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
		Read(std::string filename);
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
		/*!Stream operator that reads 2D arrays of strings*/
		Read& operator>>(Array2D<std::string>& arr);

		/*!To be used with a default constructor : opens a file named "filename", reads from the filename the type of file*/
		void open(std::string filename);

		/*!Returns the header contained in the file*/
		std::string get_header() const;
		/*!True if the file can be read*/
		bool get_status() const;

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
		bool test_ext(std::string f);

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
		size_t reading_point; //!< last bit read
};

template<typename T>
Read& Read::operator>>(T& t){
	if(unlocked){
		if(binary){ reading_point = fread(&t,sizeof(t),1,bfile); }
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
	reading_point = fread(&N,sizeof(N),1,bfile);
	if(N != mat.size()) {
		Matrice<M> mat_tmp(N);
		mat = mat_tmp;
	} 
	reading_point = fread(mat.ptr(),sizeof(M),N*N,bfile);
}

template<typename A>
void Read::read_binary_array2d(Array2D<A>& arr){
	unsigned int N_row(0);
	unsigned int N_col(0);
	reading_point = fread(&N_row,sizeof(N_row),1,bfile);
	reading_point = fread(&N_col,sizeof(N_col),1,bfile);
	if(N_row != arr.row() || N_col != arr.col()) {
		Array2D<A> arr_tmp(N_row,N_col);
		arr = arr_tmp;
	} 
	reading_point = fread(arr.ptr(),sizeof(A),N_row*N_col,bfile);
}
#endif

