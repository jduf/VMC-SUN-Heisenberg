#ifndef DEF_TRIANGLEALTERNATINGPLAQUETTE
#define DEF_TRIANGLEALTERNATINGPLAQUETTE

#include "Triangle.hpp"

class TriangleAlternatingPlaquette: public Triangle<double>{
	public:
		TriangleAlternatingPlaquette(System const& s, double const& t);
		~TriangleAlternatingPlaquette() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const t_; //!< hopping term

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;

		std::string get_mc_run_command() const;
};
#endif
