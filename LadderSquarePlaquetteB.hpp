#ifndef DEF_LADDERSQUAREPLAQUETTEB
#define DEF_LADDERSQUAREPLAQUETTEB

#include "Ladder.hpp"

class LadderSquarePlaquetteB: public Ladder<double>{
	public:
		LadderSquarePlaquetteB(System const&s, Vector<double> const& t);
		~LadderSquarePlaquetteB() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
