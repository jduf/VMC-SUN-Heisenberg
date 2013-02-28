#ifndef DEF_RST
#define DEF_RST

class Write;

#include <string>
#include <vector>
#include <cstdlib>

class RST{
	public:
		RST();
		RST(std::string filename);
		~RST();

	void title(std::string t);
		void subtitle(std::string t);
		void text(std::string t);
		void textit(std::string t);
		void textbf(std::string t);
		void item(std::string t);
		void hyperlink(std::string t, std::string l);
		void np();

		std::string RST_nl;
		std::string RST_np;
		std::string RST_title;
		std::string RST_subtitle;
		std::string RST_item;

		inline std::string get() const { return rst;};
		void set(std::string const& s) { rst = s; };

	private:
		std::string rst;
		std::vector<std::string> links;
		std::string filename;
		Write *w;
};

std::ostream& operator<<(std::ostream& flux, RST const& rst);
#endif
