#ifndef DEF_RSTfile
#define DEF_RSTfile

#include"RST.hpp"
#include"Write.hpp"

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
