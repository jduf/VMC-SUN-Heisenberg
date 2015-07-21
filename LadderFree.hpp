#ifndef DEF_LADDERFREECOMPLEX
#define DEF_LADDERFREECOMPLEX

#include "Ladder.hpp"

class LadderFree: public Ladder<double>{
	public:
		LadderFree(System const&s, Vector<double> const& t);
		~LadderFree() = default;

		void create();
		void save() const;
		void check();

	protected:
		void compute_H();
		void lattice();
		Vector<double> const t_;
};
#endif
