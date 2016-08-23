#ifndef DEF_KAGOMEPLAQUETTE6B
#define DEF_KAGOMEPLAQUETTE6B

#include "Kagome.hpp"

class KagomePlaquette6B: public Kagome<double>{
	public:
		KagomePlaquette6B(System const& s, double const& td);
		~KagomePlaquette6B() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const td_; //!< hopping terms

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
