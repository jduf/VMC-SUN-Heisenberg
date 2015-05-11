#ifndef DEF_LADDERFREECOMPLEX
#define DEF_LADDERFREECOMPLEX

#include "Ladder.hpp"

class LadderFree: public Ladder<double>{
	public:
		LadderFree(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Vector<double> const& t);
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
