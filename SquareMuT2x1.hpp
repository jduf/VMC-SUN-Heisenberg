#ifndef DEF_SQUAREMUT2X1
#define DEF_SQUAREMUT2X1

#include "Square.hpp"

class SquareMuT2x1: public Square<double>{
	public:
		SquareMuT2x1(System const& s, double const& mu, Vector<double> const& t);
		~SquareMuT2x1() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const mu_;//!< chemical potential
		Vector<double> const t_; //!< free hopping term

		void compute_H(unsigned int const& c);
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
