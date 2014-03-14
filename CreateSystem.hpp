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
		CreateSystem(unsigned int N, unsigned int n, unsigned int m, int bc, double param, Vector<unsigned int> ref);
		virtual ~CreateSystem();

		void check();
		void study(double E, double DeltaE, Vector<double> corr, std::string save_in);
		void save(Write& w) const;
		bool use_complex() const;
		bool is_bosonic() const;

		std::string get_filename() const;
		double get_param() const {return param_;}
		unsigned int get_status() const {return status_;}
		unsigned int get_N() const {return N_;}
		unsigned int get_n() const {return n_;}
		unsigned int get_m() const {return m_;}
		unsigned int get_bc() const {return bc_;}
		unsigned int get_num_links() const;
		template<typename Type>
			Matrix<Type> get_EVec() const;
		Matrix<unsigned int> get_links() const;

	private:
		unsigned int status_;
		unsigned int const N_;
		unsigned int const n_;
		unsigned int const m_;
		int const bc_;
		double param_;
		Vector<unsigned int> ref_;
		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void parse(Parseur& P);
		void create();
};
#endif
