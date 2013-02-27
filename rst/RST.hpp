#ifndef DEF_RST
#define DEF_RST

#include "Write.hpp"
#include <string>
#include <vector>
#include <cstdlib>

class RST{
	public:
		RST(std::string fname);
		~RST();

		void title(std::string t);
		void subtitle(std::string t);
		void text(std::string t);
		void item(std::string t);
		void hyperlink(std::string t, std::string l);
		void np();

		std::string RST_nl;
		std::string RST_np;
		std::string RST_title;
		std::string RST_subtitle;
		std::string RST_item;

	private:
		Write w;
		std::string fname;
		std::vector<std::string> links;
};

#endif
