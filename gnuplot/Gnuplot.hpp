#ifndef GNUPLOT
#define GNUPLOT

#include"IOFiles.hpp"
#include"Linux.hpp"

class Gnuplot {
	public:
		/*!Constructor for the creation of a .gp file in path/filename*/
		Gnuplot(std::string const&  path, std::string const& filename);
		/*!Destructor*/
		~Gnuplot(){}

		void operator=(std::string const& s){ plot_ = s + "\n"; }
		void operator+=(std::string const& s){ plot_ += s + "\n"; }

		void save_file();
		void create_image(bool silent);

		void xrange(std::string const& a, std::string const& b){plot_+="set xrange ["+a+":"+b+"]\n";}
		void xrange(double const& a, std::string const& b){xrange(tostring(a),b);}
		void xrange(std::string const& a, double const& b){xrange(a,tostring(b));}
		void xrange(double const& a, double const& b){xrange(tostring(a),tostring(b));}

		void yrange(std::string const& a, std::string const& b){plot_+="set yrange ["+a+":"+b+"]\n";}
		void yrange(double const& a, double const& b){yrange(tostring(a),tostring(b));}
		void yrange(double const& a, std::string const& b){yrange(tostring(a),b);}
		void yrange(std::string const& a, double const& b){yrange(a,tostring(b));}

		void zrange(std::string const& a, std::string const& b){plot_+="set zrange ["+a+":"+b+"]\n";}
		void zrange(double const& a, double const& b){zrange(tostring(a),tostring(b));}
		void zrange(double const& a, std::string const& b){zrange(tostring(a),b);}
		void zrange(std::string const& a, double const& b){zrange(a,tostring(b));}

	private:
		/*!Forbids copy*/
		Gnuplot(Gnuplot const& g);
		/*!Forbids assignment*/
		Gnuplot& operator=(Gnuplot g);
		
		std::string path_;		//!< path of the .gp, .png and .pdf files
		std::string filename_;	//!< filename (without the extension)
		std::string plot_;		//!< text of the .gp file
};
#endif
