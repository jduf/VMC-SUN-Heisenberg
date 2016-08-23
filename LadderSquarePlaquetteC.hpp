#ifndef DEF_LADDERSQUAREPLAQUETTEC
#define DEF_LADDERSQUAREPLAQUETTEC

#include "Ladder.hpp"

class LadderSquarePlaquetteC: public Ladder<double>{
	public:
		LadderSquarePlaquetteC(System const&s, Vector<double> const& t);
		~LadderSquarePlaquetteC() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
