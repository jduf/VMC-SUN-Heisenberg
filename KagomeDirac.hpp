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

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
