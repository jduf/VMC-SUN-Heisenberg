#ifndef DEF_LADDERDIMERA
#define DEF_LADDERDIMERA

#include "Ladder.hpp"

class LadderDimerA: public Ladder<double>{
	public:
		LadderDimerA(System const&s, Vector<double> const& t);
		~LadderDimerA() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
};
#endif
