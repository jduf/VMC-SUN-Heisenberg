#ifndef DEF_BISYSTEM
#define DEF_BISYSTEM

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include <omp.h>

class BiSystem {
	public:
		/*!Constructor that only sets param_*/
		BiSystem(Parseur& P);
		/*!Constructor that reads from file*/
		BiSystem(IOFiles& r);
		/*!Default destructor*/
		virtual ~BiSystem() = default;
		/*{Forbidden*/
		BiSystem() = delete;
		BiSystem(BiSystem const&) = delete;
		BiSystem(BiSystem&&) = delete;
		BiSystem& operator&=(BiSystem) = delete;
		/*}*/

		void add_new_param(Vector<double> const& param);
		void run(unsigned int const& nruns, unsigned int const& tmax);
		void compute_E();
		void compute_dE();

		Vector<double> const& get_E() const { return E_; }
		Matrix<double> const& get_H() const { return H_; }
		Matrix<double> const& get_O() const { return O_; }

		Vector<double> const& get_dE() const { return dE_; }
		Matrix<double> const& get_dH() const { return dH_; }
		Matrix<double> const& get_dO() const { return dO_; }

		void study();

		void save() const;

	private:
		System s_;
		std::vector<Vector<double> > param_;
		std::vector<std::vector<std::unique_ptr<MCSystem> > > mcsys_;
		Vector<double> E_;
		Matrix<double> H_;
		Matrix<double> O_;
		Vector<double> dE_;
		Matrix<double> dH_;
		Matrix<double> dO_;
};
#endif
