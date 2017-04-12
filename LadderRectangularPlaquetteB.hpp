#ifndef DEF_LADDERRECTANGULARPLAQUETTEB
#define DEF_LADDERRECTANGULARPLAQUETTEB

#include "Ladder.hpp"

class LadderRectangularPlaquetteB: public Ladder<double>{
	public:
		LadderRectangularPlaquetteB(System const&s, Vector<double> const& t);
		~LadderRectangularPlaquetteB() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
