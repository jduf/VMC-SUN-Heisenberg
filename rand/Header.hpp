#ifndef DEF_HEADER
#define DEF_HEADER

#include "Time.hpp"
#include "RST.hpp"
#include <complex>

class Header:public RST{
	public:
		Header(){};
		~Header(){};

		void init(std::string const& s);

		void add(std::string const& s);
		void add(std::string const& s, double const& d);
		void add(std::string const& s, bool const& d);
		void add(std::string const& s, unsigned int const& d);
		void add(std::string const& s, int const& d);
		void add(std::string const& s, std::string const& d);
		void add(std::string const& s, std::complex<double> const& d);
		template<typename Type>
			void add(std::string const& s, Type const& t);

	private:
		std::string when();
};

std::ostream& operator<<(std::ostream& flux, Header const& h);

template<typename Type>
void Header::add(std::string const& s, Type const& t){
	t.header_rst(s,(*this));
}
#endif
