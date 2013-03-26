#ifndef DEF_HEADER
#define DEF_HEADER

#include "Matrice.hpp"
#include "Array2D.hpp"
#include "RST.hpp"

#include <string>
#include <sstream> //-> tostring(T const& t)
#include <iostream>
#include <complex>
#include <ctime> // -> time(0)

class Header{
	public:
		Header();
		~Header();

		void init(std::string title);

		void add(std::string const& s);
		void add(std::string const& s, double const& d);
		void add(std::string const& s, bool const& d);
		void add(std::string const& s, unsigned int const& d);
		void add(std::string const& s, std::string const& d);
		void add(std::string const& s, std::complex<double> const& d);
		void add(std::string const& s, Matrice<double> const& mat);
		void add(std::string const& s, Matrice<std::complex<double> > const& mat);
		void add(std::string const& s, Array2D<unsigned int> const& arr);

		inline RST* get() const {return rst;};
		void set(std::string const& s) { rst->set(s); };
		void hyperlink(std::string const& display, std::string const& link);

	private:
		RST *rst;

		std::string when();

		template<typename T>
			std::string tostring(T const& t);
};

std::ostream& operator<<(std::ostream& flux, Header const& h);

template<typename T>
std::string Header::tostring(T const& t){
	std::ostringstream s;
	s<<t;
	return s.str();
}
#endif
