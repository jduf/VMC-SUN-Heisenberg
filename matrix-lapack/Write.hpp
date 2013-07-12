#ifndef DEF_WRITE
#define DEF_WRITE

#include "Header.hpp"

#include <fstream>

//{Description
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
//}
class Write{
	public:
		/*!Default constructor that needs a call of Write::open(std::string
		 * filename, bool binary)*/
		Write();
		/*!Opens a file named "filename", by default open a binary file*/
		Write(std::string filename);
		/*!Closes the file*/
		~Write();

		/*!To be used with a default constructor : opens a file named
		 * "filename", reads from the filename the type of file*/
		void open(std::string filename);

		/*!Stream operator that writes datas without formatting \warning
		 * doesn't work if Type=std::string => template specialization in the
		 * cpp file*/
		template<typename Type>
			Write& operator<<(Type const& t);	
		Write& operator<<(std::string const& s);
		/*!Stream operator that writes matrices, uses Matrix<Type>::operator<<
		 * \warning doesn't work if Type=std::string => template specialization
		 * in the cpp file*/
		template<typename Type>
			Write& operator<<(Matrix<Type> const& mat);
		Write& operator<<(Matrix<std::string> const& mat);
		/*!Allow to write in the header what will be written in the file*/
		template<typename Type>
			void operator()(std::string const& var, Type const& val);

		/*!Returns the filename in which the class in writing*/
		std::string get_filename() const { return filename;};

		void set_header(std::string s);

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
			void write_binary_matrix(Matrix<M> const& mat);

		std::string filename; //!< name of the file to write in
		FILE *bfile; //!< pointer on the binery file to write in
		std::ofstream tfile; //!< stream of the text file to write in
		Header *h; //!< pointer to a header (actually it will be a footer)
		bool unlocked; //!< true if the file is ready to be write in
		bool binary; //!< true if the file is binary
};

template<typename Type>//!\warning doesn't work if Type=std::string (to be checked) => template specialization in the cpp file
Write& Write::operator<<(Type const& t){
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

template<typename Type>
Write& Write::operator<<(Matrix<Type> const& mat){
	if(unlocked){
		if(binary) { write_binary_matrix(mat); }
		else { tfile<<mat<<std::flush; }
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

template<typename Type>
void Write::write_binary_matrix(Matrix<Type> const& mat){
	unsigned int N_row(mat.row());
	unsigned int N_col(mat.col());
	fwrite(&N_row, sizeof(N_row), 1 ,bfile);
	fwrite(&N_col, sizeof(N_col), 1 ,bfile);
	fwrite(mat.ptr(),sizeof(Type),mat.total(),bfile);
	fflush(bfile);
}

template<typename Type>
void Write::operator()(std::string const& var, Type const& val){
	if(h){
		(*this)<<val;
		h->add(var,val);
	} else {
		std::cerr<<"Write : operator() : there's no header in "<<filename<<std::endl;
	}	
}
#endif

