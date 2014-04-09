#ifndef DEF_RSTFILE
#define DEF_RSTFILE

#include"Write.hpp"
#include"Linux.hpp"

class RSTFile:public RST{
	public:
		RSTFile(std::string path, std::string filename);
		~RSTFile();

		void pdf() { pdf_= true; }

	private:
		std::string path_;
		std::string filename_;
		bool pdf_;
};
#endif
