#ifndef DEF_KAGOMEFERMI
#define DEF_KAGOMEFERMI

#include "Kagome.hpp"

class KagomeFermi: public Kagome<double>{
	public:
		KagomeFermi(System const& s);
		~KagomeFermi() = default;

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
