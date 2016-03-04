#ifndef DEF_KAGOMEPLAQUETTE
#define DEF_KAGOMEPLAQUETTE

#include "Kagome.hpp"

class KagomePlaquette: public Kagome<double>{
	public:
		KagomePlaquette(System const& s, double const& t);
		~KagomePlaquette() = default;

		void create();
		void check();

	protected:
		double const t_;

		void compute_H();
		void display_results();
		void lattice();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
