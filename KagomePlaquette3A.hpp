#ifndef DEF_KAGOMEPLAQUETTE3A
#define DEF_KAGOMEPLAQUETTE3A

#include "Kagome.hpp"

class KagomePlaquette3A: public Kagome<double>{
	public:
		KagomePlaquette3A(System const& s, double const& td);
		~KagomePlaquette3A() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		double const td_;

		void compute_H();
		void display_results();
		void lattice();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
