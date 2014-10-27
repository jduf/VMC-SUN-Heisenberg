#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "ChainFermi.hpp"
#include "ChainPolymerized.hpp"

#include "TriangleFermi.hpp"

#include "SquareFermi.hpp"
#include "SquarePiFlux.hpp"
#include "SquareJastrow.hpp"

#include "KagomeFermi.hpp"
#include "KagomeDirac.hpp"
#include "KagomeVBC.hpp"

#include "Honeycomb0pp.hpp"

#include "Parseur.hpp"
#include "List.hpp"

class CreateSystem{
	public:
		CreateSystem(Parseur& P);
		CreateSystem(IOFiles* r);
		virtual ~CreateSystem();

		void init(IOFiles* read=NULL, IOSystem* ios=NULL);
		bool is_over() const { return over_; }

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
		void create();
		/*}*/

		/*{Other class calls*/
		/*!Calls Fermionic::is_degenerate()*/
		bool is_degenerate() const {
			if(RGL_){return RGL_->is_degenerate();}
			if(CGL_){return CGL_->is_degenerate();}
			return true;
		}
		/*{!Calls System::get_status()
		 *status_ = 10: GenericSystem<Type> is not defined
		 * status_ = 3 : System is initialized (in System)
		 * status_ = 2 : System is allowed (in System)
		 * status_ = 1 : System*D is allowed (in System*D)
		 * status_ = 0 : Found an initial state (in SystemFermionic)
		 */
		/*}*/
		unsigned int get_status() const {
			if(RGL_){return RGL_->get_status();}
			if(CGL_){return CGL_->get_status();}
			return 10;
		}
		/*}*/

		/*!Returns a pointer on the GenericSystem created*/
		System const* get_system() const { 
			if(RGL_){return RGL_;}
			if(CGL_){return CGL_;}
			return NULL;
		}
		/*!Returns false if ref_(1)==1 : real eigenvectors*/
		bool use_complex() const {
			if(ref_(1) == 1){ return false; }
			else { return true; }
		}
		/*!Returns true if ref_(1)==0 : Jastrow wavefunction*/
		bool is_bosonic() const {
			if(ref_(1) == 0){ return true; }
			else { return false; }
		}

	private:
		/*!Forbids copy*/
		CreateSystem(CreateSystem const& cs);
		/*!Forbids assignment*/
		CreateSystem& operator=(CreateSystem cs);

		Vector<unsigned int> ref_;
		unsigned int const N_;
		unsigned int const m_;
		unsigned int const n_;
		Vector<unsigned int> M_;
		int const bc_;
		List<double> d_;
		unsigned int type_;
		bool over_;
		unsigned int sel0_;
		unsigned int sel1_;

		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void parse(Parseur& P);
		void error();
};
#endif
