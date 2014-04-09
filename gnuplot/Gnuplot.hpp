#ifndef GNUPLOT
#define GNUPLOT

#include"Write.hpp"
#include"Linux.hpp"

class Gnuplot {
	public:
		Gnuplot(std::string path, std::string filename, std::string type);
		Gnuplot(std::string path, std::string filename);
		~Gnuplot();

		void save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Matrix<double> const& z);
		void save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y, Vector<double> const& z);
		void save_data(std::string data_file, Vector<double> const& x, Vector<double> const& y);
		void save_data(std::string data_file, Matrix<double> const& z);
		void preplot(std::string s);
		void add(std::string s);

		void save_file();
		void create_image(bool silent);

	private:
		std::string path_;
		std::string filename_;
		std::string plottype_;
		std::string preplot_;
		std::string plot_;
};
#endif
