#ifndef DEF_HEADER
#define DEF_HEADER

#include "Read.hpp"
#include "Write.hpp"
#include "Matrice.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <complex>
#include <ctime>

struct RSTFormatting{
	RSTFormatting():item("+ "),bf("**"),it("*"),np("\n\n"),title('='){}
	std::string item,bf,it,np;
	char title;
};

class Header{
	public:
		Header(Read& r);
		Header(std::string const& s);

		void add(std::string str, double const& d);
		void add(std::string str, std::complex<double> const& c);

		void show(){std::cout<<h;}
		void write(Write& w);
		void add(std::string str, Matrice<double>& m);

		inline std::string get() const{ return h;};

	private:
		std::string h;
		RSTFormatting rst;

		std::string title(std::string& t);

		void when();
		void intro(std::string const& s);

		template<typename T>
			std::string tostring(T const& t);
};

template<typename T>
std::string Header::tostring(T const& t){
	std::ostringstream s;
	s<<t;
	return s.str();
}

#endif
