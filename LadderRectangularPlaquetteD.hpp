#ifndef DEF_LADDERRECTANGULARPLAQUETTED
#define DEF_LADDERRECTANGULARPLAQUETTED

#include "Ladder.hpp"

class LadderRectangularPlaquetteD: public Ladder<double>{
	public:
		LadderRectangularPlaquetteD(System const&s, Vector<double> const& t);
		~LadderRectangularPlaquetteD() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
