#ifndef DEF_IOFILES
#define DEF_IOFILES

#include "Header.hpp"
#include <fstream>

/*{Description*/
/*!Class that allows to read/write datas from/in a binary or a text file.
 *
 *When used with a text file
 * - the datas must be separated with a space
 * - only one matrix or array is safely read per file
 * - the size of the matrix or array has to be specified
 *
 *When used with a binary_ file
 * - many different kinds of datas can be extracted
*/
/*}*/
class IOFiles{
	public:
		/*!Opens a file named "filename_", reads from the filename_ the type of file*/
		IOFiles(std::string filename, bool write);
		/*!Closes the file*/
		~IOFiles();

		/*!Stream operator forall types 
		 * \warning doesn't work if Type=std::string */
		template<typename Type>
			IOFiles& operator>>(MyType<Type>& t);
		IOFiles& operator>>(double& t){read_basic_type(t); return (*this);}
		IOFiles& operator>>(std::complex<double>& t){read_basic_type(t); return (*this);}
		IOFiles& operator>>(int& t){read_basic_type(t); return (*this);}
		IOFiles& operator>>(unsigned int& t){read_basic_type(t); return (*this);}
		IOFiles& operator>>(bool & t){read_basic_type(t); return (*this);}
		IOFiles& operator>>(std::string& t){read_string(t); return (*this);}

		template<typename Type>
			IOFiles& operator<<(MyType<Type> const& t);
		IOFiles& operator<<(double const& t){write_basic_type(t); return (*this);}
		IOFiles& operator<<(std::complex<double> const& t){write_basic_type(t); return (*this);}
		IOFiles& operator<<(int const& t){write_basic_type(t); return (*this);}
		IOFiles& operator<<(unsigned int const& t){write_basic_type(t); return (*this);}
		IOFiles& operator<<(bool const& t){write_basic_type(t); return (*this);}
		IOFiles& operator<<(std::string const& t){
			write_string(t); return (*this);}

		template<typename Type>
			void read(Type* m, unsigned int const& N);
		template<typename Type>
			void read(Type& m);

		template<typename Type>
			void write(Type* m, unsigned int const& N);
		template<typename Type>
			void write(Type const& m);

		template<typename Type>
			void operator()(std::string const& var, Type const& val);

		/*!Returns true if the file is open as a binary file*/
		bool is_binary() const {return binary_;}
		/*!Returns true if the file is open*/
		bool is_open() const {return open_;}
		/*!Returns the filename_ in which the class in writing*/
		std::string get_filename() const { return filename_;};

		/*!Returns the header contained in the file*/
		std::string get_header() const;
		void add_to_header(std::string const& s);

		static std::string endl; //!<Gives a way to end lines

	private:
		/*!Forbids copy constructor*/
		IOFiles(IOFiles const& R);
		/*!Forbids assertion operator*/
		IOFiles& operator=(IOFiles const&);

		/*!Subroutine that check if the correct extension is given*/
		void test_ext();

		/*!Subroutine needed to open a binary file*/
		void open_binary();
		/*!Subroutine needed to open a text file*/
		void open_txt();
		/*!Extract the header*/
		void read_header();
		/*!Write the header*/
		void write_header();

		template<typename Type>
			void read_basic_type(Type& t);
		template<typename Type>
			void write_basic_type(Type const& t);

		void read_string(std::string& t);
		void write_string(std::string const& t);

		std::string filename_;	//!< name of the file to read from
		bool write_;
		bool binary_; 			//!< true if the file is binary_
		bool open_;				//!< true if the file is ready to be read from
		std::fstream file_;		//!< text file to read form
		Header* h_;				//!< pointer to a header (actually it will be a footer)
};

template<typename Type>
IOFiles& IOFiles::operator>>(MyType<Type>& t){
	if(open_ && !write_){
		if(binary_){ t.read_from_file(*this);}
		else{file_>>t;}
	} else {
		std::cerr<<"IOFiles::operator>>(MyType<Type>) : can't read from "<<filename_<<std::endl;
	}
	return (*this);
}

template<typename Type>
IOFiles& IOFiles::operator<<(MyType<Type> const& t){
	if(open_ && write_){
		if(binary_){t.write_in_file(*this);}
		else{file_<<t<<std::flush;}
	} else {
		std::cerr<<"IOFiles::operator<<(MyType<Type>) : can't write in "<<filename_<<std::endl;
	}
	return (*this);
}

template<typename Type>
void IOFiles::read_basic_type(Type& t){
	if(open_ && !write_){
		if (binary_){read(t);}
		else {file_>>t;}
	} else {
		std::cerr<<"IOFiles::read_basic_type(Type) : can't read from "<<filename_<<std::endl;
	}
}

template<typename Type>
void IOFiles::write_basic_type(Type const& t){
	if(open_ && write_){
		if (binary_){write(t);}
		else {file_<<t;}
	} else {
		std::cerr<<"IOFiles::write_basic_type(Type) : can't write in "<<filename_<<std::endl;
	}
}

template<typename Type>
void IOFiles::read(Type* m, unsigned int const& N){
	file_.read((char*)(m),N*sizeof(Type));
}

template<typename Type>
void IOFiles::write(Type* m, unsigned int const& N){
	file_.write((char*)(m),N*sizeof(Type));
}

template<typename Type>
void IOFiles::read(Type& m){
	file_.read((char*)(&m),sizeof(Type));
}

template<typename Type>
void IOFiles::write(Type const& m){
	file_.write((char *)(&m),sizeof(Type));
}

template<typename Type>
void IOFiles::operator()(std::string const& var, Type const& val){
	if(h_ && open_ && write_){
		(*this)<<val;
		h_->add(var,val);
	} else {
		std::cerr<<"IOFiles::operator()(string,val) : can't write in "<<filename_<<std::endl;
	}	
}
#endif


