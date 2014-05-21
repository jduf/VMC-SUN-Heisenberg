#ifndef GNUPLOT
#define GNUPLOT

#include"IOFiles.hpp"
#include"Linux.hpp"

class Gnuplot {
	public:
		Gnuplot(std::string const&  path, std::string const& filename);
		~Gnuplot();

		void operator=(std::string const& s){ plot_ = s + "\n"; }
		void operator+=(std::string const& s){ plot_ += s + "\n"; }

		void save_file();
		void create_image(bool silent);

		void xrange(std::string const& a, std::string const& b){plot_+="set xrange ["+a+":"+b+"]\n";}
		void xrange(double const& a, double const& b){xrange(tostring(a),tostring(b));}
		void xrange(double const& a, std::string const& b){xrange(tostring(a),b);}
		void xrange(std::string const& a, double const& b){xrange(a,tostring(b));}

		void yrange(std::string const& a, std::string const& b){plot_+="set yrange ["+a+":"+b+"]\n";}
		void yrange(double const& a, double const& b){yrange(tostring(a),tostring(b));}
		void yrange(double const& a, std::string const& b){yrange(tostring(a),b);}
		void yrange(std::string const& a, double const& b){yrange(a,tostring(b));}

	private:
		std::string path_;
		std::string filename_;
		std::string plot_;
};
#endif
