#ifndef DEF_CHAINPOLYMERIZEDJJP
#define DEF_CHAINPOLYMERIZEDJJP

#include "ChainPolymerized.hpp"

class ChainPolymerizedJJp: public ChainPolymerized {
	public:
		ChainPolymerizedJJp(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, Vector<double> const& t, Vector<double> const& JJp);
		~ChainPolymerizedJJp(){}
		
		void save() const;

	private:
		Vector<double> JJp_;//!< polymerization parameter 

		std::string extract_level_7();
		std::string extract_level_6();
};
#endif

