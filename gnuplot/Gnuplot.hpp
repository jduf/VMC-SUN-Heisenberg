#ifndef GNUPLOT
#define GNUPLOT

#include"Matrix.hpp"
#include"Write.hpp"

class Gnuplot {
	public:
		Gnuplot(std::string filename, std::string type);
		~Gnuplot();

		void data(Matrix<double> const& x, Matrix<double> const& y, Matrix<double> const& z);
		void code(std::string s);
		void save();

		void test();

	private:
		Matrix<double> gp_data;
		std::string filename;
		std::string gp_code;
};

#endif
