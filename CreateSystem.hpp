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
#include "SquareACSL.hpp"
#include "SquareFreeComplex.hpp"

#include "KagomeFermi.hpp"
#include "KagomeDirac.hpp"
#include "KagomeVBC.hpp"

#include "Honeycomb0pp.hpp"

class CreateSystem{
	public:
		/*{Description*/
		/*!Takes a pointer to an already existing System (better than to create
		 * a System within this class because in VMCMinimization, many
		 * GenericSystem will be created therefore if one can avoid the
		 * re-creation of a System, it saves ressources) */
		/*}*/
		CreateSystem(System const* const s);
		virtual ~CreateSystem();
		/*{Forbidden*/
		CreateSystem() = delete;
		CreateSystem(CreateSystem const&) = delete;
		CreateSystem(CreateSystem&&) = delete;
		CreateSystem& operator=(CreateSystem cs) = delete;
		/*}*/

		/*{Description*/
		/*!If read!=NULL, it will use parameters saved in IOFiles instead of
		 * those stored in Container C_*/
		/*}*/
		void init(Vector<double> const* const param, Container* C);
		void create(bool const& try_solve_degeneracy=false);

		/*{IOSystem calls*/
		/*!The attributes of IOSystem will be copied to RGL_/CGL_*/
		void set_IOSystem(IOSystem const* const ios){
			if(RGL_){ RGL_->set_IOSystem(ios); }
			if(CGL_){ CGL_->set_IOSystem(ios); }
		}
		/*!Calls IOSystem::analyse(unsigned int const& level)*/
		std::string analyse(unsigned int const& level) {
			if(RGL_){ return RGL_->analyse(level); }
			if(CGL_){ return CGL_->analyse(level); }
			return "";
		}
		/*!Returns the filename (only usefull for mc) */
		std::string get_filename() const {
			if(RGL_){ return RGL_->get_filename(); }
			if(CGL_){ return CGL_->get_filename(); }
			return "";
		}
		/*!Returns the path (only usefull for mc) */
		std::string get_path() const {
			if(RGL_){ return RGL_->get_path();}
			if(CGL_){ return CGL_->get_path();}
			return "";
		}
		/*}*/

		/*{GenericSystem calls*/
		/*!Calls GenericSystem::save_input() virtual method*/
		void save_param(IOFiles& w) const {
			if(RGL_){ RGL_->save_param(w); }
			if(CGL_){ CGL_->save_param(w); }
		}
		/*!Calls GenericSystem::check() pure virtual method*/
		void check() const {
			if(RGL_){ return RGL_->check(); }
			if(CGL_){ return CGL_->check(); }
		}
		/*!Calls GenericSystem::create() pure virtual method*/
		void get_wf_symmetries(std::vector<Matrix<int> >& sym) const {
			if(RGL_){ RGL_->get_wf_symmetries(sym); }
			if(CGL_){ CGL_->get_wf_symmetries(sym); }
		}
		void lattice(std::string const& path, std::string const& filename) const {
			if(RGL_){ RGL_->lattice(path,filename); }
			if(CGL_){ CGL_->lattice(path,filename); }
		}
		/*}*/

		/*{System calls*/
		/*!Calls System::get_status() : see System.hpp*/
		unsigned int get_status() const {
			if(RGL_){ return RGL_->get_status(); }
			if(CGL_){ return CGL_->get_status(); }
			return 10;
		}
		/*!Calls System::get_status() : see System.hpp*/
		void set_bonds(System* const s) const {
			if(RGL_){ s->set_bonds(RGL_); }
			if(CGL_){ s->set_bonds(CGL_); }
		}
		/*}*/

		/*!Returns ref*/
		Vector<unsigned int> const&  get_ref() const { return ref_; }
		/*!Returns a pointer on the GenericSystem created*/
		System const* get_GS() const {
			if(RGL_){ return RGL_; }
			if(CGL_){ return CGL_; }
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
		Container C_;
		System const* s_;
		Vector<unsigned int> ref_;

		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void error() const;
};
#endif
