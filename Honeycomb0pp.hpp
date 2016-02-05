#ifndef DEF_HONEYCOMB0PP
#define DEF_HONEYCOMB0PP

#include "Honeycomb.hpp"

class Honeycomb0pp: public Honeycomb<double>{
	public:
		Honeycomb0pp(System const& s, double const& td);
		~Honeycomb0pp() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		double const td_;

		void compute_H();

		void display_results();
		void lattice();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;

		std::string extract_level_7();
};
#endif
