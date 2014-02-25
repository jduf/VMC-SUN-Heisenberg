#ifndef GNUPLOT
#define GNUPLOT

#include"Write.hpp"

class Gnuplot {
	public:
		Gnuplot(std::string filename, std::string type);
		~Gnuplot();

		void save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Matrix<double> const& z);
		void save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Vector<double> const& z);
		void save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y);
		void save_data(std::string data_file, Matrix<double> const& z);
		void add_plot_param(std::string s);
		void preplot(std::string s);

		void test();

	private:
		std::string filename_;
		std::string plottype_;
		std::string preplot_;
		std::string plot_;
};
#endif
