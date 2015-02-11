#ifndef DEF_RST
#define DEF_RST

#include <string>
#include <iostream>

class RST{
	public:
		/*!Constructor*/
		RST();
		/*!Constructor that sets rst_=rst*/
		RST(std::string rst);
		/*!Destructor*/
		virtual ~RST(){}

		void title(std::string t, std::string symb);
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

		std::string const& get() const { return rst_;};
		void set(std::string const& s="") { rst_ = s; };

	protected:
		std::string RST_nl_;	//!< string for a new line
		std::string RST_np_;	//!< string for a new paragraph
		std::string RST_item_;	//!< string for a new item
		std::string rst_;		//!< text of the .rst file
};

std::ostream& operator<<(std::ostream& flux, RST const& rst);

#include <sstream> //-> tostring(T const& t)
template<typename Type>
std::string tostring(Type const& t){
	std::ostringstream s;
	s<<t;
	return s.str();
}
#endif
