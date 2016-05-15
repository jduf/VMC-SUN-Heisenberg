#ifndef DEF_KAGOMEPLAQUETTE3B
#define DEF_KAGOMEPLAQUETTE3B

#include "Kagome.hpp"

class KagomePlaquette3B: public Kagome<double>{
	public:
		KagomePlaquette3B(System const& s, double const& td);
		~KagomePlaquette3B() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const td_;

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
