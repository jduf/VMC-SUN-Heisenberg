#ifndef DEF_LADDERFREECOMPLEX
#define DEF_LADDERFREECOMPLEX

#include "Ladder.hpp"

class LadderFree: public Ladder<double>{
	public:
		LadderFree(System const&s, Vector<double> const& t);
		~LadderFree() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();
		void get_wf_symmetries(std::vector<Matrix<int> >& sym) const;

	protected:
		Vector<double> const t_;

		void compute_H();
		void lattice();

		std::string extract_level_6();
		void plot();

	private:
		unsigned int set_spuc(Vector<double> const& t);
};
#endif
