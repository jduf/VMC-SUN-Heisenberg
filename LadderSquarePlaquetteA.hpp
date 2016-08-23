#ifndef DEF_LADDERSQUAREPLAQUETTEA
#define DEF_LADDERSQUAREPLAQUETTEA

#include "Ladder.hpp"

class LadderSquarePlaquetteA: public Ladder<double>{
	public:
		LadderSquarePlaquetteA(System const&s, Vector<double> const& t);
		~LadderSquarePlaquetteA() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
