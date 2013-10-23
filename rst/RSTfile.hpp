#ifndef DEF_RSTFILE
#define DEF_RSTFILE

#include"Write.hpp"
#include"Linux.hpp"

class RSTfile:public RST{
	public:
		RSTfile(std::string filename, std::string path="");
		~RSTfile();

		void pdf(bool c_pdf) { create_pdf_=c_pdf; }

	private:
		std::string path_;
		std::string filename_;
		bool create_pdf_;
};
#endif
