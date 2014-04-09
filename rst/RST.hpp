#ifndef DEF_RST
#define DEF_RST

#include <string>
#include <sstream> //-> tostring(T const& t)
#include <iostream>

class RST{
	public:
		RST();
		RST(std::string rst);
		virtual ~RST();

		void title(std::string t,std::string symb);
		void text(std::string t);
		void textit(std::string t);
		void textbf(std::string t);
		void item(std::string t);
		void lineblock(std::string t);
		void def(std::string t, std::string def);
		void hyperlink(std::string display, std::string link);
		void figure(std::string image, std::string legend, unsigned int width=1000);
		void link_figure(std::string image, std::string legend, std::string link, unsigned int width=1000);
		void np();
		void nl();

		inline std::string get() const { return rst;};
		inline void set(std::string const& s="") { rst = s; };

	protected:
		std::string RST_nl;
		std::string RST_np;
		std::string RST_item;
		std::string rst;
};

template<typename Type>
std::string tostring(Type const& t){
	std::ostringstream s;
	s<<t;
	return s.str();
}

std::ostream& operator<<(std::ostream& flux, RST const& rst);
#endif
