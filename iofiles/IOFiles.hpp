#ifndef DEF_IOFILES
#define DEF_IOFILES

#include "Header.hpp"
#include <fstream>
#include <cstring>

/*{*/
 /*!Class that allows to read/write datas from/in a binary or a text file.
 * 
 * When used with a text file
 *  - the datas must be separated with a space
 *  - only one matrix or array is safely read per file
 *  - the size of the matrix or array has to be specified
 * 
 * When used with a binary_ file
 *  - many different kinds of datas can be extracted */
/*}*/
class IOFiles{
	public:
		/*{Description*/
		/*!Opens a file named filename_ and deduces from filename_ if the file
		 * is a binary or a text file*/
		/*}*/
		IOFiles(std::string filename, bool write);
		/*!Destructor that closes the file*/
		~IOFiles();

		/*!Stream operator for all types 
		 * \warning doesn't work if Type=std::string */
		IOFiles& operator>>(double& t){read(&t,1,sizeof(double)); return (*this);}
		IOFiles& operator>>(std::complex<double>& t){read(&t,1,sizeof(std::complex<double>)); return (*this);}
		IOFiles& operator>>(int& t){read(&t,1,sizeof(int)); return (*this);}
		IOFiles& operator>>(unsigned int& t){read(&t,1,sizeof(unsigned int)); return (*this);}
		IOFiles& operator>>(bool & t){read(&t,1,sizeof(bool)); return (*this);}
		IOFiles& operator>>(std::string& t){read_string(t); return (*this);}

		IOFiles& operator<<(double const& t){write(&t,1,sizeof(double)); return (*this);}
		IOFiles& operator<<(std::complex<double> const& t){write(&t,1,sizeof(std::complex<double>)); return (*this);}
		IOFiles& operator<<(int const& t){write(&t,1,sizeof(int)); return (*this);}
		IOFiles& operator<<(unsigned int const& t){write(&t,1,sizeof(unsigned int)); return (*this);}
		IOFiles& operator<<(bool const& t){write(&t,1,sizeof(bool)); return (*this);}
		IOFiles& operator<<(const char* t){write_string(t,std::strlen(t)); return (*this);}
		IOFiles& operator<<(std::string const& t){write_string(t.c_str(),t.size()); return (*this);}

		/*Allow to extract personal classes or data*/
		template<typename Type>
			void read(Type* m, unsigned int const& N, size_t const& type_size);
		/*!Read file and extract a Type value. Useful for initialization*/
		template<typename Type>
			Type read();
		/*Allow to write personal classes or data*/
		template<typename Type>
			void write(Type* m, unsigned int const& N, size_t const& type_size);
		/*Write val with label var*/
		template<typename Type>
			void write(std::string const& var, Type const& val);

		/*!Returns file_*/
		std::fstream& stream(){ return file_;}
		/*!Returns true if the file is open as a binary file*/
		bool is_binary() const {return binary_;}
		/*!Returns true if the file is open*/
		bool is_open() const {return open_;}
		/*!Returns the filename_ in which the class in writing*/
		std::string get_filename() const { return filename_;};

		/*!Change the precision on the output text files*/
		void precision(unsigned int const& N);
		/*!Returns the header contained in the file*/
		std::string get_header() const;
		/*!Returns the header contained in the file*/
		Header* add_header() const { return header_; };

		static std::string const endl; //!<Gives a way to end lines

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

		/*!Read string*/
		void read_string(std::string& t);
		/*!Write string*/
		void write_string(const char* t, unsigned int const& N);

		std::string filename_;	//!< name of the file to read from
		bool write_;			//!< true if the file is writable
		bool binary_; 			//!< true if the file is binary_
		bool open_;				//!< true if the file is ready to be read from
		std::fstream file_;		//!< text file to read form
		Header* header_;		//!< pointer to a header (actually it will be a footer)
};

template<typename Type>
void IOFiles::read(Type* t, unsigned int const& N, size_t const& type_size){
	if(open_ && !write_){
		if (binary_){ file_.read((char*)(t),N*type_size); }
		else { file_>>*t; }
	} else {
		std::cerr<<"IOFiles::read(Type*,unsigned int,size_t) : can't read from "<<filename_<<std::endl;
	}
}

template<typename Type>
void IOFiles::write(Type* t, unsigned int const& N, size_t const& type_size){
	if(open_ && write_){
		if (binary_){ file_.write((char*)(t),N*type_size); }
		else { file_<<*t; }
	} else {
		std::cerr<<"IOFiles::write(Type*,unsigned int,size_t) : can't write in "<<filename_<<std::endl;
	}
}

template<typename Type>
Type IOFiles::read(){
	Type t;
	(*this)>>t;
	return t;
}

template<typename Type>
void IOFiles::write(std::string const& var, Type const& val){
	if(header_ && open_ && write_){
		(*this)<<val;
		header_->add(var,val);
	} else {
		std::cerr<<"void IOFiles::write(string,val) : can't write in "<<filename_<<std::endl;
	}	
}
#endif
