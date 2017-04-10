#ifndef DEF_LADDERRECTANGULARPLAQUETTEC
#define DEF_LADDERRECTANGULARPLAQUETTEC

#include "Ladder.hpp"

class LadderRectangularPlaquetteC: public Ladder<double>{
	public:
		LadderRectangularPlaquetteC(System const&s, Vector<double> const& t);
		~LadderRectangularPlaquetteC() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		unsigned int set_spuc();

		void display_results();
};
#endif
