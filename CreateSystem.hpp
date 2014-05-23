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

		void create(double const& x, unsigned int const& type){
			if(RGL_){RGL_->create(x,type);}
			if(CGL_){CGL_->create(x,type);}
		}

		System const& get_system() const { 
			if(RGL_) {return RGL_->get_system();}
			else {return CGL_->get_system();}
		}
		std::string get_filename() const {
			if(RGL_) {return RGL_->get_filename();}
			else { return CGL_->get_filename();}
		}
		bool use_complex() const {
			if(ref_(1) == 1){ return false; }
			else { return true; }
		}
		bool is_bosonic() const {
			if(ref_(1) == 0){ return true; }
			else { return false; }
		}
		bool is_degenerate() const {
			if(RGL_) {return RGL_->is_degenerate();}
			if(CGL_) {return CGL_->is_degenerate();}
			return true;
		}
		template<typename Type>
			Bosonic<Type> const& get_bosonic() const;
		template<typename Type>
			Fermionic<Type> const& get_fermionic() const;

		void save(IOFiles& w) const{
			if(RGL_){RGL_->save(w);}
			if(CGL_){CGL_->save(w);}
		}

		void check(){
			if(RGL_){return RGL_->check();}
			if(CGL_){return CGL_->check();}
		}

	private:
		Vector<unsigned int> ref_;
		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void parse(Parseur& P);
		void init(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc);
};
#endif
