#ifndef DEF_MCSIM
#define DEF_MCSIM

#include "CreateSystem.hpp"
#include "MonteCarlo.hpp"

class MCSim {
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
		void copy_S(std::unique_ptr<MCSystem> const& MCS);
		/*!Creates MonteCarlo, then run on MCS_*/
		void run(unsigned int const& thermalization_steps, unsigned int const& tmax);
		/*}*/

		/*{System and MCSystem calls*/
		/*!Calls void System::set_observable(unsigned int const& which)*/
		void set_observables(unsigned int const& which){ MCS_->set_observables(which); }
		/*!Calls bool System::check_conv(double const& convergence_criterion)*/
		bool check_conv(double const& convergence_criterion){
			return MCS_->check_conv(convergence_criterion);
		}
		/*!Calls void System::complete_analysis(double const& convergence_criterion)*/
		void complete_analysis(double const& convergence_criterion){
			MCS_->complete_analysis(convergence_criterion);
		}
		/*!Calls virtual void MCSystem::free_memory() = 0*/
		void free_memory(){ MCS_->free_memory(); }
		/*}*/

		/*{Write in IOFiles methods*/
		void write(IOFiles& w) const;
		/*!Save the result in a single file (wavefunction parameters, observables)*/
		void save(IOFiles& w) const;
		/*}*/

		/*{Static methods*/
		static bool sort_by_E(MCSim const& a, MCSim const& b);
		static unsigned int sort_by_param_for_merge(MCSim const& list, MCSim const& new_elem);
		static void merge(MCSim& list, MCSim& new_elem);
		/*}*/

		/*{Simple value return*/
		/*!Returns param_*/
		Vector<double> const& get_param() const { return param_; }
		/*!Returns MCS_*/
		std::unique_ptr<MCSystem> const& get_MCS() const { return MCS_; }
		/*!Returns true if MCS_ can be run by MonteCarlo*/
		bool is_created() const { return (MCS_.get() && !MCS_->get_status()); }
		/*}*/

	private:
		Vector<double> param_;
		std::unique_ptr<MCSystem> MCS_;
};
#endif
