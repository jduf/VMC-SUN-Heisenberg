#ifndef DEF_RSTFILE
#define DEF_RSTFILE

#include"IOFiles.hpp"
#include"Linux.hpp"

class RSTFile:public RST{
	public:
		/*!Constructor for the creation of a .rst file in path/filename*/
		RSTFile(std::string const& path, std::string const& filename);
		/*!Destructor*/
		~RSTFile(){}

		/*!Saves the .rst file and create a .html and a .pdf if pdf==true*/
		void save(bool pdf);

	private:
		std::string path_;		//!< path of the .rst, .html and .pdf files
		std::string filename_;	//!< filename (without the extension)
};
#endif
