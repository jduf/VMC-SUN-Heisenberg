#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "ChainFermi.hpp"
#include "ChainFree.hpp"
#include "ChainPolymerized.hpp"

#include "LadderFermi.hpp"
#include "LadderFree.hpp"
#include "LadderFreeFlux.hpp"

#include "TriangleFermi.hpp"
#include "TriangleFree.hpp"
#include "TrianglePlaquette.hpp"
#include "TriangleMu.hpp"
#include "TrianglePhi.hpp"
#include "TriangleChiral.hpp"

#include "SquareFermi.hpp"
#include "SquareFree.hpp"
#include "SquareMu.hpp"
#include "SquareDimerizedBar.hpp"
#include "SquareFreeFlux.hpp"
#include "SquarePiFlux.hpp"
#include "SquareChiral.hpp"
#include "SquareJastrow.hpp"

#include "KagomeFermi.hpp"
#include "KagomeDirac.hpp"
#include "KagomeVBC.hpp"

#include "HoneycombFermi.hpp"
#include "HoneycombFree.hpp"
#include "Honeycomb0pp.hpp"
#include "HoneycombPiFlux.hpp"
#include "HoneycombChiral.hpp"

/*{*//*!Class that creates any kind of wavefunctions
	   and gives a way to act on it (set observables, merge simulations,
	   compute errors on measurments, save...) 

	   This class handles everything that is required to create a wavefunction
	   from a child of GenergicSystem. 

	   It can also store (merge) and handle (complete_analysis, save,
	   display_results,..) the results of the VMC simulation contained in a
	   child of MCSystem.*//*}*/
class CreateSystem{
	public:
		/*{*//*!Takes a pointer to an already existing System 
			   (better than to create a System within this class because in
			   VMCMinimization, many GenericSystem will be created therefore if
			   one can avoid the re-creation of a System, it saves ressources)
			   *//*}*/
		CreateSystem(System const* const s);
		virtual ~CreateSystem();
		/*{Forbidden*/
		CreateSystem() = delete;
		CreateSystem(CreateSystem const&) = delete;
		CreateSystem(CreateSystem&&) = delete;
		CreateSystem& operator=(CreateSystem) = delete;
		/*}*/

		/*{Core methods*/
		void init(Vector<double> const* const param, Container* C);
		void create(bool const& try_solve_degeneracy=false);
		/*}*/

		/*{System calls*/
		/*!Calls System::merge : see System*/
		void merge(System* const s) const {
			if(RGL_){ return RGL_->merge(s); }
			if(CGL_){ return CGL_->merge(s); }
		}
		/*!Calls System::complete_analysis : see System*/
		void complete_analysis(double const& convergence_criterion) const {
			if(RGL_){ RGL_->complete_analysis(convergence_criterion); }
			if(CGL_){ CGL_->complete_analysis(convergence_criterion); }
		}
		/*}*/

		/*{IOSystem calls*/
		/*!The attributes of IOSystem will be copied to RGL_/CGL_*/
		void set_IOSystem(IOSystem const* const ios){
			if(RGL_){ RGL_->set_IOSystem(ios); }
			if(CGL_){ CGL_->set_IOSystem(ios); }
		}
		/*!Calls IOSystem::get_system_info*/
		RST const& get_system_info() const {
			if(RGL_){ return RGL_->get_system_info(); }
			else    { return CGL_->get_system_info(); }
		}
		/*!Calls IOSystem::analyse*/
		std::string analyse(unsigned int const& level){
			if(RGL_){ return RGL_->analyse(level); }
			else    { return CGL_->analyse(level); }
		}
		/*}*/

		/*{GenericSystem calls*/
		/*!Calls GenericSystem::set_obs pure virtual method*/
		void set_obs(int const& nobs) const {
			if(RGL_){ return RGL_->set_obs(nobs); }
			if(CGL_){ return CGL_->set_obs(nobs); }
		}
		/*!Calls GenericSystem::get_wf_symmetries pure virtual method*/
		void get_wf_symmetries(std::vector<Matrix<int> >& sym) const {
			if(RGL_){ RGL_->get_wf_symmetries(sym); }
			if(CGL_){ CGL_->get_wf_symmetries(sym); }
		}
		/*!Calls GenericSystem::check pure virtual method*/
		void check() const {
			if(RGL_){ return RGL_->check(); }
			if(CGL_){ return CGL_->check(); }
		}
		/*!Calls GenericSystem::display_results pure virtual method*/
		void display_results() const {
			if(RGL_){ RGL_->display_results(); }
			if(CGL_){ CGL_->display_results(); }
		}
		/*}*/

		/*{Simple value return*/
		/*!Returns a pointer on the GenericSystem created*/
		System const* get_GenericSystem() const {
			if(RGL_){ return RGL_; }
			if(CGL_){ return CGL_; }
			return NULL;
		}
		/*!Returns true if ref_(1)==2 : complex eigenvectors*/
		bool use_complex() const { return ref_(1)==2; }
		/*!Returns true if ref_(1)==0 : Jastrow wavefunction*/
		bool is_bosonic() const { return ref_(1)==0; }
		/*!Returns the filename*/
		std::string const& get_filename() const {
			if(RGL_){ return RGL_->get_filename(); }
			else    { return CGL_->get_filename(); }
		}
		/*!Returns the path*/
		std::string const& get_path() const {
			if(RGL_){ return RGL_->get_path(); }
			else    { return CGL_->get_path(); }
		}
		/*!Calls System::get_status : see System*/
		unsigned int get_status() const {
			if(RGL_){ return RGL_->get_status(); }
			if(CGL_){ return CGL_->get_status(); }
			return 10;
		}
		/*!Calls System::get_obs : see System*/
		std::vector<Observable> const& get_obs() const {
			if(RGL_){ return RGL_->get_obs(); }
			else { return CGL_->get_obs(); }
		}
		/*}*/

		/*{Print and output in IOFiles methods*/
		/*!Calls GenericSystem::save_param and System::save*/
		void save(IOFiles& w) const {
			if(RGL_){
				RGL_->save_param(w); 
				RGL_->save(w); 
			}
			if(CGL_){ 
				CGL_->save_param(w); 
				CGL_->save(w); 
			}
		}
		/*!Calls System::print*/
		void print(unsigned int const& nobs) const {
			if(RGL_){ RGL_->print(nobs); }
			if(CGL_){ CGL_->print(nobs); }
		}
		/*}*/

	private:
		Container* C_;
		System const* s_;
		Vector<unsigned int> ref_;

		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		/*Returns an error message when a wavefunction isn't found*/
		void error() const {
			std::cerr<<__PRETTY_FUNCTION__<<" : ref_ = ["<<ref_<<"] unknown"<<std::endl;
		}
};
#endif
