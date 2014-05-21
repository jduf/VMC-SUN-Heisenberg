#ifndef DEF_PSTRICKS
#define DEF_PSTRICKS

#include "Linux.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

class PSTricks {
	public:
		PSTricks(std::string path, std::string filename, bool silent);
		~PSTricks();

		void add(std::string s);
		void line(std::string linetype, double x0, double y0, double x1, double y1, std::string options);
		void frame(double x0, double y0, double x1, double y1, std::string options="");
		void put(double x, double y, std::string s, std::string options="");
		void pie(Vector<double> const& x, double r, std::string options="");
		void polygon(Matrix<double> const& x, std::string options="");

		void pdf(){ pdf_ = true;}

	private:
		std::string path_;
		std::string filename_;
		std::string s_;
		bool pdf_;
		bool silent_;
};
#endif
