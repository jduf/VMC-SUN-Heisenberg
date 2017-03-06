#ifndef DEF_SQUAREK2MU2X2
#define DEF_SQUAREK2MU2X2

#include "Square.hpp"

class Squarek2Mu2x2: public Square<double>{
	public:
		Squarek2Mu2x2(System const& s, double const& mu, Vector<double> const& t);
		~Squarek2Mu2x2() = default;

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
