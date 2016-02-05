#ifndef DEF_HONEYCOMBFERMI
#define DEF_HONEYCOMBFERMI

#include "Honeycomb.hpp"

class HoneycombFermi: public Honeycomb<double>{
	public:
		HoneycombFermi(System const& s);
		~HoneycombFermi() = default;

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
