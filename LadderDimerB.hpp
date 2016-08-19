#ifndef DEF_LADDERDIMERB
#define DEF_LADDERDIMERB

#include "Ladder.hpp"

class LadderDimerB: public Ladder<double>{
	public:
		LadderDimerB(System const&s, Vector<double> const& t);
		~LadderDimerB() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
