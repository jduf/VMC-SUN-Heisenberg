#ifndef DEF_MCSIM
#define DEF_MCSIM

#include "CreateSystem.hpp"
#include "MonteCarlo.hpp"

class MCSim{
	public:
		/*!Constructor that only sets param_*/
		MCSim(Vector<double> const& param);
		/*!Constructor that reads from file*/
		MCSim(IOFiles& r);
		/*!Default destructor*/
		virtual ~MCSim() = default;
		/*{Forbidden*/
		MCSim() = delete;
		MCSim(MCSim const&) = delete;
		MCSim(MCSim&&) = delete;
		MCSim& operator&=(MCSim) = delete;
		/*}*/

		/*{Core methods*/
		/*!Sets MCS_ to a new MCSystem created via C*/
		void create_S(System const* const s);
		/*!Sets MCS_ to a copy obtained via MCSystem::clone() run on MCS*/
		void copy_S(std::shared_ptr<MCSim> const& mcsim);
		/*!Creates MonteCarlo, then run on MCS_*/
		void run(unsigned int const& ts, unsigned int const& tmax);
		/*}*/

		/*{System and MCSystem calls*/
		/*!Calls void System::set_obs(Observable const& obs)*/
		void set_obs(Observable const& obs){ MCS_->set_obs(obs); }
		/*!Calls void System::set_obs(unsigned int const& i)*/
		void set_obs(unsigned int const& i){ MCS_->set_obs(i); }

		/*!Calls bool System::check_conv(double const& convergence_criterion)*/
		bool check_conv(double const& convergence_criterion){ return MCS_->check_conv(convergence_criterion); }
		/*!Calls void System::complete_analysis(double const& convergence_criterion)*/
		void complete_analysis(double const& convergence_criterion){ MCS_->complete_analysis(convergence_criterion); }
		/*!Calls void System::merge(System* const s)*/
		void merge(std::shared_ptr<MCSim> const& mcsim){ MCS_->merge(mcsim->MCS_.get()); }

		/*!Calls virtual void MCSystem::free_memory() = 0*/
		void free_memory(){ MCS_->free_memory(); }
		/*!Calls System::print()*/
		void print(unsigned int const& nobs) const { MCS_->print(nobs); }
		/*}*/

		/*{Output in IOFiles methods*/
		/*!Write raw data (no output in header) => made to same many MCSim*/
		void write(IOFiles& w) const;
		/*!Write nice data (with output header) => made to same one MCSim*/
		void save(IOFiles& w) const;
		/*}*/

		/*{Static methods*/
		static bool sort_by_E(MCSim const& a, MCSim const& b);
		static unsigned int sort_for_merge(MCSim const& list, MCSim const& new_elem);
		static unsigned int sort_by_param_for_merge(Vector<double> const& a, Vector<double> const& b);
		/*}*/

		/*{Simple value return*/
		/*!Return the Data<double> Energy*/
		Data<double> const& get_energy() const { return MCS_->get_energy(); }
		/*!Returns param_*/
		Vector<double> const& get_param() const { return param_; }
		/*!Returns true if MCS_ can be run by MonteCarlo*/
		bool is_created() const { return (MCS_.get() && !MCS_->get_status()); }
		/*}*/

	private:
		Vector<double> param_;			//!< variational parameters used to create the wavefunction
		std::unique_ptr<MCSystem> MCS_; //!< pointer on a MCSystem that can be measured via MonteCarlo
};
#endif
