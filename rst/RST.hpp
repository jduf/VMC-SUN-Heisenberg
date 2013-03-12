#ifndef DEF_RST
#define DEF_RST

class Write;

#include <string>
//#include <vector>
#include <cstdlib> // -> system(std::string command)

class RST{
	public:
		RST();
		RST(std::string filename);
		~RST();

		void title(std::string t,std::string symb);
		void text(std::string t);
		void textit(std::string t);
		void textbf(std::string t);
		void item(std::string t);
		void def(std::string t, std::string def);
		void hyperlink(std::string display, std::string link);
		void np();

		std::string RST_nl;
		std::string RST_np;
		std::string RST_item;

		inline std::string get() const { return rst;};
		void set(std::string const& s) { rst = s; };

	private:
		std::string rst;
		//std::vector<std::string> links;
		std::string filename;
		Write *w;
};

std::ostream& operator<<(std::ostream& flux, RST const& rst);
#endif
