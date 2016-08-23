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
		IOFiles(std::string const& filename, bool const& write, bool const& append_txt_file);
		/*!Destructor that closes the file*/
		~IOFiles();
		/*{Forbidden*/
		IOFiles() = delete;
		IOFiles(IOFiles const&) = delete;
		IOFiles(IOFiles&&) = delete;
		IOFiles& operator=(IOFiles&) = delete;
		/*}*/

		/*!Stream operator for all types
		 * \warning needs an overload for std::string*/
		template<typename Type>
			IOFiles& operator>>(Type& t){ read(&t,1,sizeof(Type)); return (*this); }
		IOFiles& operator>>(std::string& t);
		/*Allows to extract personal classes or data*/
		template<typename Type>
			void read(Type* m, unsigned int const& N, size_t const& type_size);
		/*!Reads file and extract a Type value. Useful for initialization*/
		template<typename Type>
			Type read();

		/*!Stream operator for all types
		 * \warning needs an overload for std::string*/
		template<typename Type>
			IOFiles& operator<<(Type const& t){ write(&t,1,sizeof(Type)); return (*this); }
		IOFiles& operator<<(std::string const& t);
		/*Allows to write personal classes or data*/
		template<typename Type>
			void write(Type* m, unsigned int const& N, size_t const& type_size);
		/*Writes val with label var*/
		template<typename Type>
			void write(std::string const& var, Type const& val);

		/*!Returns file_*/
		std::fstream& stream(){ return file_; }
		/*!Returns true if the file is open as a binary file*/
		bool const& is_binary() const { return binary_; }
		/*!Returns true if the file is open*/
		bool const& is_open() const { return open_; }
		/*!Returns the filename_ in which the class in writing*/
		std::string const& get_filename() const { return filename_; }

		/*!Changes the precision on the output text files*/
		void precision(unsigned int const& N);
		/*!Returns the header contained in the file*/
		std::string get_header() const;
		/*!Returns the header contained in the file*/
		Header* add_header() const { return header_; }

		static std::string const endl; //!<Gives a way to end lines

	private:
		/*!Extracts the header*/
		void read_header();
		/*!Writes the header*/
		void write_header();

		std::string filename_;	//!< name of the file to read from
		bool write_;			//!< true if the file is writable
		bool binary_; 			//!< true if the file is binary_
		Header* header_;		//!< pointer to a header (actually it will be a footer)
		std::fstream file_;		//!< text file to read form
		bool open_;				//!< true if the file is ready to be read from
};

template<typename Type>
void IOFiles::read(Type* t, unsigned int const& N, size_t const& type_size){
	if(open_ && !write_){
		if (binary_){ file_.read((char*)(t),N*type_size); }
		else { file_>>*t; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : can't read from "<<filename_<<std::endl; }
}

template<typename Type>
Type IOFiles::read(){
	Type t;
	(*this)>>t;
	return t;
}

template<typename Type>
void IOFiles::write(Type* t, unsigned int const& N, size_t const& type_size){
	if(open_ && write_){
		if (binary_){ file_.write((char*)(t),N*type_size); }
		else { file_<<*t; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : can't write in "<<filename_<<std::endl; }
}

template<typename Type>
void IOFiles::write(std::string const& var, Type const& val){
	if(header_ && open_ && write_){
		(*this)<<val;
		header_->add(var,val);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : can't write in "<<filename_<<std::endl; }	
}
#endif
