#ifndef DEF_PSTRICKS
#define DEF_PSTRICKS

#include "Write.hpp"
#include "Linux.hpp"

class PSTricks {
	public:
		PSTricks(std::string filename, std::string path="");
		~PSTricks();

		void add(std::string s);
		void line(std::string linetype, double x0, double y0, double x1, double y1, std::string options="");
		void frame(double x0, double y0, double x1, double y1, std::string options="");
		void put(double x, double y, std::string s);

	private:
		std::string path_;
		std::string filename_;
		std::string s_;

};
#endif
