#ifndef DEF_KAGOMEVBC
#define DEF_KAGOMEVBC

#include "Kagome.hpp"

class KagomeVBC: public Kagome<std::complex<double> >{
	public:
		KagomeVBC(System const& s);
		~KagomeVBC() = default;

		void create();
		void check();

	private:
		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
