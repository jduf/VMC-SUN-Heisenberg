#ifndef DEF_RST
#define DEF_RST

#include <string>
#include <cstdlib> // -> system(std::string command)

class RST{
	public:
		RST();
		~RST();

		void title(std::string t,std::string symb);
		void text(std::string t);
		void textit(std::string t);
		void textbf(std::string t);
		void item(std::string t);
		void lineblock(std::string t);
		void def(std::string t, std::string def);
		void hyperlink(std::string display, std::string link);
		void figure(std::string image, std::string legend, unsigned int scale=100);
		void np();

		std::string RST_nl;
		std::string RST_np;
		std::string RST_item;

		inline std::string get() const { return rst;};
		inline void set(std::string const& s) { rst = s; };

	protected:
		std::string rst;
};

std::ostream& operator<<(std::ostream& flux, RST const& rst);
#endif
