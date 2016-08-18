#ifndef DEF_PSTRICKS
#define DEF_PSTRICKS

#include "Linux.hpp"
#include "Vector.hpp"

class PSTricks {
	public:
		PSTricks(std::string path, std::string filename);
		/*!Default destructor*/
		~PSTricks() = default;
		/*{Forbidden*/
		PSTricks() = delete;
		PSTricks(PSTricks const&) = delete;
		PSTricks(PSTricks&&) = delete;
		PSTricks operator=(PSTricks) = delete;
		/*}*/

		void begin(double const& xbl, double const& ybl, double const& xtr, double const& ytr, std::string const& imagename);
		void end(bool const& silent, bool const& png, bool const& crop);

		void add(std::string const& s);
		void line(std::string const& linetype, double const& x0, double const& y0, double const& x1, double const& y1, std::string const& options);
		void linked_lines(std::string const& linetype, Matrix<double> const& x, std::string const& options="");
		void frame(double const& x0, double const& y0, double const& x1, double const& y1, std::string const& options="");
		void put(double const& x, double const& y, std::string const& s, std::string const& options="");
		void pie(double const& x, double const& y, Vector<double> const& p, double const& r, std::string const& options="");
		void circle(Vector<double> const& x, double const& r, std::string const& options="");
		void cross(Vector<double> const& x, double const& r, std::string const& options="");
		void polygon(Matrix<double> const& x, std::string const& options="");

	private:
		std::string path_;
		std::string filename_;
		std::string s_;
		bool pdf_;
		bool silent_;
		bool begin_end_;
};
#endif
