#ifndef DEF_RST
#define DEF_RST

#include "Miscellaneous.hpp"

class RST{
	public:
		/*!Default constructor*/
		RST() = default;
		/*!Constructor that sets rst_=rst*/
		RST(std::string const& rst);
		/*!Default destructor*/
		virtual ~RST(){};
		/*{Forbidden*/
		RST(RST const&) = delete;
		RST(RST&&) = delete;
		RST& operator=(RST) = delete;
		/*}*/

		void title(std::string const& t, std::string const& symb);
		void text(std::string const& t);
		void math_centered(std::string const& t);
		void item(std::string const& t);
		void lineblock(std::string const& t);
		void def(std::string const& t, std::string const& def);
		void hyperlink(std::string const& display, std::string const& link);
		void figure(std::string const& image, std::string const& legend, unsigned int width=1000);
		void link_figure(std::string const& image, std::string const& legend, std::string const& link, unsigned int width=1000);
		void np();
		void nl();

		static std::string textit(std::string const& t);
		static std::string textbf(std::string const& t);
		static std::string math(std::string const& t);

		std::string const& get() const { return rst_;};
		void set(std::string const& s="") { rst_ = s; };

		static std::string const nl_;	//!< string for a new line
		static std::string const np_;	//!< string for a new paragraph
		static std::string const item_;	//!< string for a new item

	protected:
		std::string rst_;//!< text of the .rst file
};

std::ostream& operator<<(std::ostream& flux, RST const& rst);
#endif
