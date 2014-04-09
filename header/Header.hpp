#ifndef DEF_HEADER
#define DEF_HEADER

#include "Matrix.hpp"
#include "Vector.hpp"
#include "Time.hpp"
#include "RST.hpp"

class Header{
	public:
		Header();
		~Header();

		void init(std::string title);

		void add(std::string const& s);
		void add(std::string const& s, double const& d);
		void add(std::string const& s, bool const& d);
		void add(std::string const& s, unsigned int const& d);
		void add(std::string const& s, int const& d);
		void add(std::string const& s, std::string const& d);
		void add(std::string const& s, std::complex<double> const& d);
		void add(std::string const& s, Vector<unsigned int> const& vec);
		void add(std::string const& s, Vector<double> const& vec);
		void add(std::string const& s, Matrix<int> const& mat);
		void add(std::string const& s, Matrix<unsigned int> const& mat);
		void add(std::string const& s, Matrix<double> const& mat);
		void add(std::string const& s, Matrix<std::complex<double> > const& mat);

		void set(std::string const& s) { rst.set(s); };
		RST get() const {return rst;};

	private:
		RST rst;

		std::string when();
};

std::ostream& operator<<(std::ostream& flux, Header const& h);
#endif
