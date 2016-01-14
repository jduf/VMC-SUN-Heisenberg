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

		Matrix<double> const& get_H() const { return H_; }
		Matrix<double> const& get_dH() const { return dH_; }
		Matrix<double> const& get_O() const { return O_; }
		Matrix<double> const& get_dO() const { return dO_; }

		Vector<double> const& get_E() const { return E_; }
		Vector<double> const& get_dE() const { return dE_; }

		void save(IOFiles& w) const;

	private:
		System s_;
		std::vector<Vector<double> > param_;
		std::vector<std::vector<std::unique_ptr<MCSystem> > > mcsys_;
		Matrix<double> H_;
		Matrix<double> O_;
		Matrix<double> dH_;
		Matrix<double> dO_;
		Vector<double> E_;
		Vector<double> dE_;
};
#endif
