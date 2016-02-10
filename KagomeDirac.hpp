#ifndef DEF_KAGOMEDIRAC
#define DEF_KAGOMEDIRAC

#include "Kagome.hpp"

class KagomeDirac: public Kagome<double>{
	public:
		KagomeDirac(System const& s);
		~KagomeDirac() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();
		void lattice();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
