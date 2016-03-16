#ifndef DEF_KAGOMECHIRAL
#define DEF_KAGOMECHIRAL

#include "Kagome.hpp"

class KagomeChiral: public Kagome<std::complex<double> >{
	public:
		KagomeChiral(System const& s, double const& phi);
		~KagomeChiral() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		double const phi_; //!< flux per triangular plaquette

		void compute_H();
		void display_results();
		void lattice();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
