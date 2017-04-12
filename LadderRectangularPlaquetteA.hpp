#ifndef DEF_LADDERRECTANGULARPLAQUETTEA
#define DEF_LADDERRECTANGULARPLAQUETTEA

#include "Ladder.hpp"

class LadderRectangularPlaquetteA: public Ladder<double>{
	public:
		LadderRectangularPlaquetteA(System const&s, Vector<double> const& t);
		~LadderRectangularPlaquetteA() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
