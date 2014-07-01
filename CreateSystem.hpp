#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "ChainPolymerized.hpp"
#include "ChainFermi.hpp"
#include "SquarePiFlux.hpp"
#include "KagomeFermi.hpp"
#include "KagomeDirac.hpp"
#include "KagomeVBC.hpp"
#include "Parseur.hpp"
#include "List.hpp"

class CreateSystem{
	public:
		CreateSystem(Parseur& P);
		virtual ~CreateSystem();

		void init();
		void create(){
			if(RGL_){RGL_->create();}
			if(CGL_){CGL_->create();}
		}

		System const* get_system() const { 
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

		void save(IOFiles& w) const{
			if(RGL_){RGL_->save(w);}
			if(CGL_){CGL_->save(w);}
		}

		void check(){
			if(RGL_){return RGL_->check();}
			if(CGL_){return CGL_->check();}
		}

		bool is_over(){ return over_; }

	private:
		CreateSystem(CreateSystem const& cs);

		Vector<unsigned int> ref_;
		unsigned int const N_;
		unsigned int const n_;
		unsigned int const m_;
		Vector<unsigned int> M_;
		int const bc_;
		List<double> d_;
		unsigned int type_;
		bool over_;

		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void parse(Parseur& P);
};
#endif
