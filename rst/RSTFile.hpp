#ifndef DEF_RSTFILE
#define DEF_RSTFILE

#include"IOFiles.hpp"
#include"Linux.hpp"

class RSTFile:public RST{
	public:
		RSTFile(std::string path, std::string filename);
		~RSTFile();

		void save(bool pdf);

	private:
		std::string path_;
		std::string filename_;
};
#endif
