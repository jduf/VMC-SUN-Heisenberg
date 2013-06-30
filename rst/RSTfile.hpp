#ifndef DEF_RSTFILE
#define DEF_RSTFILE

#include"RST.hpp"
#include"Write.hpp"
#include"Linux.hpp"

class RSTfile:public RST{
	public:
		RSTfile(std::string filename, std::string path="");
		~RSTfile();

		void pdf(bool c_pdf) { create_pdf=c_pdf; }

	private:
		std::string path;
		std::string filename;
		Write w;
		bool create_pdf;
};
#endif
