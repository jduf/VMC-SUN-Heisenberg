#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "ChainFermi.hpp"
#include "ChainDimerized.hpp"
#include "ChainTrimerized.hpp"

#include "SquareFermi.hpp"
#include "SquareMu.hpp"
#include "SquarePiFlux.hpp"
#include "SquareSU2PhiFlux.hpp"
#include "SquareJastrow.hpp"

#include "HoneycombSU4.hpp"
#include "HoneycombSU3.hpp"

#include "TriangleFermi.hpp"
#include "TriangleMu.hpp"
#include "TrianglePhi.hpp"
#include "TriangleJastrow.hpp"

#include "Parseur.hpp"

class CreateSystem{
	public:
		CreateSystem(Parseur& P);
		virtual ~CreateSystem();

		void create();
		void save();

	private:
		Container param_;
		Vector<unsigned int> ref_;
		GenericSystem<double>* Sr_;
		GenericSystem<std::complex<double> >* Sc_;
		bool valid_;
};
#endif
