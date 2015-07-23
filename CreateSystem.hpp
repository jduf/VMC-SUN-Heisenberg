#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "ChainFermi.hpp"
#include "ChainPolymerized.hpp"

#include "LadderFermi.hpp"
#include "LadderFree.hpp"

#include "TriangleFermi.hpp"

#include "SquareFermi.hpp"
#include "SquarePiFlux.hpp"
#include "SquareJastrow.hpp"
//#include "SquareFreeReal.hpp"
#include "SquareACSL.hpp"
#include "SquareFreeComplex.hpp"

//#include "KagomeFermi.hpp"
//#include "KagomeDirac.hpp"
//#include "KagomeVBC.hpp"

#include "Honeycomb0pp.hpp"

#include "Parseur.hpp"

class CreateSystem{
	public:
		CreateSystem(System const* const s);
		virtual ~CreateSystem();
		/*{Forbidden*/
		CreateSystem() = delete;
		CreateSystem(CreateSystem const&) = delete;
		CreateSystem(CreateSystem&&) = delete;
		CreateSystem& operator=(CreateSystem cs) = delete;
		/*}*/

		void set_param(Container* C, Vector<double> const* param = NULL);
		void init(IOFiles* read=NULL, IOSystem* ios=NULL);

		/*{IOSystem calls*/
		/*!Calls IOSystem::analyse(unsigned int const& level)*/
		std::string analyse(unsigned int const& level) {
			if(RGL_){return RGL_->analyse(level);}
			if(CGL_){return CGL_->analyse(level);}
			return "";
		}
		/*!Calls IOSystem::init_output_file(output)*/
		void init_output_file(IOFiles& output) const {
			if(RGL_){RGL_->init_output_file(output);}
			if(CGL_){CGL_->init_output_file(output);}
		}
		/*!Returns the filename (only usefull for mc) */
		std::string get_filename() const {
			if(RGL_){return RGL_->get_filename();}
			if(CGL_){return CGL_->get_filename();}
			return "";
		}
		/*!Returns the path (only usefull for mc) */
		std::string get_path() const {
			if(RGL_){return RGL_->get_path();}
			if(CGL_){return CGL_->get_path();}
			return "";
		}
		/*}*/

		/*{GenericSystem calls*/
		/*!Calls GenericSystem::save() virtual method*/
		void save() const {
			if(RGL_){RGL_->save();}
			if(CGL_){CGL_->save();}
		}
		/*!Calls GenericSystem::check() pure virtual method*/
		void check() const {
			if(RGL_){return RGL_->check();}
			if(CGL_){return CGL_->check();}
		}
		/*!Calls GenericSystem::create() pure virtual method*/
		void create(bool try_solve_degeneracy=false);
		/*}*/

		/*{System calls*/
		/*!Calls System::get_status() : see System.hpp*/
		unsigned int get_status() const {
			if(RGL_){return RGL_->get_status();}
			if(CGL_){return CGL_->get_status();}
			return 10;
		}
		/*!Calls System::get_status() : see System.hpp*/
		void set_bonds(System* const s) const {
			if(RGL_){ s->set_bonds(RGL_); }
			if(CGL_){ s->set_bonds(CGL_); }
		}
		/*}*/

		/*!Returns ref*/
		static Vector<unsigned int> get_ref(std::string const& wf);
		Vector<unsigned int> const&  get_ref() const { return ref_; }
		/*!Returns a pointer on the GenericSystem created*/
		System const* get_system() const { 
			if(RGL_){return RGL_;}
			if(CGL_){return CGL_;}
			return NULL;
		}
		/*!Returns true if ref_(1)==2 : complex eigenvectors*/
		bool use_complex() const {
			if(ref_(1) == 2){ return true; }
			else { return false; }
		}
		/*!Returns true if ref_(1)==0 : Jastrow wavefunction*/
		bool is_bosonic() const {
			if(ref_(1) == 0){ return true; }
			else { return false; }
		}

	private:
		System const* s_;
		Container C_;
		Vector<unsigned int> ref_;

		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void error() const;
};
#endif
