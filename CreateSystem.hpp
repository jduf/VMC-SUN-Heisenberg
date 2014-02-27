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

		void create();
		void save(Write& w);
		bool use_complex() const;
		bool is_bosonic() const;

		unsigned int get_N() const {return N_;}
		unsigned int get_n() const {return n_;}
		unsigned int get_m() const {return m_;}
		unsigned int get_num_links() const;
		std::string get_filename() const;
		Matrix<unsigned int> get_sts() const;
		template<typename Type>
			Matrix<Type> get_EVec() const;

	private:
		unsigned int N_;
		unsigned int n_;
		unsigned int m_;
		Vector<unsigned int> ref_;
		GenericSystem<double>* RGL_;
		GenericSystem<std::complex<double> >* CGL_;

		void parse(Parseur& P);
};
#endif
