#ifndef DEF_TRIANGLEMU
#define DEF_TRIANGLEMU

#include "Triangle.hpp"

class TriangleMu: public Triangle<double>{
	public:
		TriangleMu(System const& s, double const& mu);
		~TriangleMu() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const mu_;//!< chemical potential

		void compute_H(unsigned int const& c);
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;

		std::string get_mc_run_command() const;
};
#endif
