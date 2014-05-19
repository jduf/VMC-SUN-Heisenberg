#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "ChainPolymerized.hpp"
#include "ChainFermi.hpp"
#include "SquarePiFlux.hpp"
#include "Parseur.hpp"

class CreateSystem{
	public:
		CreateSystem(Parseur& P);
		CreateSystem(CreateSystem const& cs, double param);
		virtual ~CreateSystem();

		void check();
		void save(IOFiles& w) const;

		std::string get_filename() const {
			if(RGL_){return RGL_->get_filename();}
			if(CGL_){return CGL_->get_filename();}
			return 0;
		}
		bool use_complex() const {
			if(ref_(1) == 1){ return false; }
			else { return true; }
		}
		bool is_bosonic() const {
			if(ref_(1) == 0){ return true; }
			else { return false; }
		}
		bool ready() const {return ready_;}
		double get_param() const {return param_;}
		unsigned int get_N() const {return N_;}
		unsigned int get_n() const {return n_;}
		unsigned int get_m() const {return m_;}
		unsigned int get_bc() const {return bc_;}
		unsigned int get_num_links() const {
			if(RGL_){return RGL_->get_num_links();}
			if(CGL_){return CGL_->get_num_links();}
			return 0;
		}
		Matrix<unsigned int> get_links() const {
			if(RGL_){return RGL_->get_links();}
			if(CGL_){return CGL_->get_links();}
			return 0;
		}
		template<typename Type>
			Matrix<Type> get_EVec() const;

	private:
		unsigned int const N_;
		unsigned int const n_;
		unsigned int const m_;
		double const param_;
		int const bc_;
		bool ready_;
		Vector<unsigned int> ref_;
		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void parse(Parseur& P);
		void create();
};
#endif
