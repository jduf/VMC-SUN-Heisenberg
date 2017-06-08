#ifndef DEF_SQUAREK2MUPI
#define DEF_SQUAREK2MUPI

#include "Square.hpp"

class Squarek2MuPi: public Square<double>{
	public:
		Squarek2MuPi(System const& s, double const& mu);
		~Squarek2MuPi() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const mu_;//!< chemical potential

		void compute_H(unsigned int const& c);
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
