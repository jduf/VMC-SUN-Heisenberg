#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "ChainPolymerized.hpp"
#include "ChainFermi.hpp"
#include "SquarePiFlux.hpp"
#include "Parseur.hpp"

class CreateSystem{
	public:
		CreateSystem(Parseur& P);
		virtual ~CreateSystem();

		void check();
		void save(IOFiles& w) const;
		void create(double const& x);

		bool use_complex() const {
		if(ref_(1) == 1){ return false; }
		else { return true; }
		}
		bool is_bosonic() const {
			if(ref_(1) == 0){ return true; }
			else { return false; }
		}

		System* get_system() const { 
			if(RGL_){return RGL_;}
			if(CGL_){return CGL_;}
			return NULL;
		}
		template<typename Type>
		Bosonic<Type>* get_bosonic() const;
		template<typename Type>
		Fermionic<Type>* get_fermionic() const;

	private:
		Vector<unsigned int> ref_;
		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void parse(Parseur& P);
		void init(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc);
};
#endif
