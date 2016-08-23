#ifndef DEF_LADDERFERMI
#define DEF_LADDERFERMI

#include "Ladder.hpp"

/*{Description*/
/*!Creates a ladder with uniform hopping parameter */
/*}*/
class LadderFermi: public Ladder<double>{
	public:
		LadderFermi(System const& s);
		~LadderFermi() = default;

		void create();
		void check();

	private:
		void compute_H();
		void display_results();
};
#endif
