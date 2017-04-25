#ifndef DEF_CHAINPOLYMERIZEDJJP
#define DEF_CHAINPOLYMERIZEDJJP

#include "Chain.hpp"

class ChainPolymerizedJJp: public Chain<double> {
	public:
		ChainPolymerizedJJp(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, Vector<double> const& t, Vector<double> const& JJp);
		~ChainPolymerizedJJp() = default;

		void create();
		void save() const;
		void check();

	private:
		Vector<double> t_;//!< polymerization parameter
		Vector<double> JJp_;//!< polymerization parameter

		void compute_H(unsigned int const& c);
		std::string extract_level_7();
		std::string extract_level_6();
};
#endif
