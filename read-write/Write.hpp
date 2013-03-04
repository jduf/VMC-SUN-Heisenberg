#ifndef DEF_WRITE
#define DEF_WRITE

#include "Matrice.hpp"
#include "Array2D.hpp"
#include "Header.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <complex>

/*!Class that allows to write datas easily in a file.
 *To be used with Read.hpp
 *Can be used to write a binary or a text file.
 *
 *When used with a text file
 * - the datas must be separated with a space
 * - only one matrix or array is safely saved per file
 *
 *When used with a binary file
 * - many different kinds of datas can be saved together
*/
class Write{
	public:
		/*!Default constructor that needs a call of Write::open(std::string filename, bool binary)*/
		Write();
		/*!Opens a file named "filename", by default open a binary file*/
		Write(std::string filename);
		/*!Closes the file*/
		~Write();

		/*!Stream operator that writes datas without formatting*/
		template<typename T>
			Write& operator<<(T const& t);	
		/*!Stream operator that writes matrices, uses Matrice<M>::operator<<*/
		template<typename M>
			Write& operator<<(Matrice<M> const& mat);
		/*!Stream operator that writes arrays, uses Array2D<A>::operator<<*/
		template<typename A>
			Write& operator<<(Array2D<A> const& arr);
		/*!Stream operator that writes strings*/
		Write& operator<<(std::string const& s);
		Write& operator<<(Array2D<std::string> const& arr);

		/*!To be used with a default constructor : opens a file named "filename", reads from the filename the type of file*/
		void open(std::string filename);

		template<typename T>
			void operator()(std::string const& var, T const& val);

		/*!Returns the filename in which the class in writing*/
		std::string get_filename() const { return filename;};

		static std::string endl; //!<Gives a way to end lines

	private:
		/*!Forbids copy constructor*/
		Write(Write const& s);
		/*!Forbids assertion operator*/
		Write& operator=(Write const&);

		/*!Subroutine needed to open a binary file*/
		void open_binary();
		/*!Subroutine needed to open a text file*/
		void open_txt();
		/*!Write the filename and the date in the header*/
		void write_header();

		/*!Subroutine that check if the correct extension is given*/
		bool test_ext(std::string f);

		/*!Subroutine needed to write a matrix in a binary file*/
		template<typename M>
			void write_binary_matrix(Matrice<M> const& mat);
		/*!Subroutine needed to write a matrix in a binary file*/
		template<typename A>
			void write_binary_array2d(Array2D<A> const& arr);

		std::string filename; //!< name of the file to write in
		FILE *bfile; //!< pointer on the binery file to write in
		std::ofstream tfile; //!< stream of the text file to write in
		Header *h; //!< pointer to a header (actually it will be a footer)
		bool unlocked; //!< true if the file is ready to be write in
		bool binary; //!< true if the file is binary
};

template<typename T>//doesn't work with string (it seems)
Write& Write::operator<<(T const& t){
	if(unlocked){
		if(binary){
			fwrite(&t, sizeof(t), 1 ,bfile);
			fflush(bfile);
		} else {
			tfile<< t<<std::flush; 
		}
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename M>//doesn't work with string
Write& Write::operator<<(Matrice<M> const& mat){
	if(unlocked){
		if(binary) { write_binary_matrix(mat); }
		else { tfile<<mat<<std::flush; }
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename A>
Write& Write::operator<<(Array2D<A> const& arr){
	if(unlocked){
		if(binary) { write_binary_array2d(arr);}
		else { tfile<<arr<<std::flush; }
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename M>//doesn't work with string
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

template<typename A>//doesn't work with string
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

template<typename T>
void Write::operator()(std::string const& var, T const& val){
	if(h){
		(*this)<< val;
		h->add(var,val);
	} else {
		std::cout<<"Write : there's no header in "<<filename<<std::endl;
	}	
}
#endif

