#ifndef GNUPLOT
#define GNUPLOT

#include"Matrix.hpp"
#include"Write.hpp"

class Gnuplot {
	public:
		Gnuplot(std::string filename, std::string type);
		~Gnuplot();

		void save_data(std::string data_file, Matrix<double> const& x, Matrix<double> const& y, Matrix<double> const& z);
		void save_data(std::string data_file, Matrix<double> const& x, Matrix<double> const& y);
		void code(std::string s);
		void save_code();

		void test();

	private:
		std::string filename;
		std::string gp_code;
};
#endif
