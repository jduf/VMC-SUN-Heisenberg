#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Parseur.hpp"

#include "ChainDimerized.hpp"
#include "ChainFermi.hpp"
#include "ChainPolymerized.hpp"
#include "ChainTrimerized.hpp"
#include "HoneycombSU3.hpp"
#include "HoneycombSU4.hpp"
#include "SquareFermi.hpp"
#include "SquareJastrow.hpp"
#include "SquareMu.hpp"
#include "SquarePiFlux.hpp"
#include "SquareSU2PhiFlux.hpp"
#include "TriangleFermi.hpp"
#include "TriangleJastrow.hpp"
#include "TriangleMu.hpp"
#include "TrianglePhi.hpp"

class CreateSystem{
	public:
		CreateSystem(Parseur& P);
		CreateSystem(Container const& C);
		virtual ~CreateSystem();

		void create();
		void save();
		void get_param(Container& param);
		void get_input(Container& input);
		bool use_complex();

	private:
		Container param_;
		GenericSystem<double>* Sr_;
		GenericSystem<std::complex<double> >* Sc_;
		Vector<unsigned int> ref_;

		void parse(Parseur& P, Container& C);
};
#endif
