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
		/*!Sets MCS_ to a new MCSystem created via a pointer to System*/
		void create_S(System const* const s);
		/*!Sets MCS_ to a copy obtained via MCSystem::clone() run on MCS*/
		void copy_clear_S(std::shared_ptr<MCSim> const& mcsim);
		/*!Creates MonteCarlo, then run on MCS_*/
		void run(unsigned int const& ts, unsigned int const& tmax);
		/*!Creates a CreateSystem and call its display_results method*/
		void display_results(std::string const& filename, std::string const& sim, std::string const& info, std::string const& analyse, std::string const& path, std::string const& dir, RSTFile* const rst_file, bool const& replace_title_with_link_in_rst);
		/*!Creates a CreateSystem and call its display_results method with default filenames and no RST*/
		void display_results();
		/*}*/

		/*{System and MCSystem calls*/
		/*!Calls void System::set_obs(Observable const& obs)*/
		void set_obs(Observable const& obs){ MCS_->set_obs(obs); }
		/*!Gets relative incertitude of the energy*/
		double get_dEoE(){ return MCS_->get_dEoE(); }
		/*!Calls bool System::check_conv(double const& convergence_criterion)*/
		bool check_conv(){ return MCS_->check_conv(); }
		/*!Calls void System::complete_analysis(double const& convergence_criterion)*/
		void complete_analysis(){ MCS_->complete_analysis(); }
		/*!Calls void System::merge(System* const s)*/
		void merge(std::shared_ptr<MCSim> const& mcsim){ MCS_->merge(mcsim->MCS_.get()); }
		/*!Calls virtual void MCSystem::free_memory() = 0*/
		void free_memory(){ MCS_->free_memory(); }
		/*!Calls System::print(bool const& all)*/
		void print(bool const& all) const { MCS_->print(all); }
		/*}*/

		/*{Output in IOFiles methods*/
		/*!Write raw data (no output in header): save many MCSim*/
		void write(IOFiles& w) const;
		/*!Write formatted data (with output header): save one MCSim or in text files*/
		void save(IOFiles& w) const;
		/*}*/

		/*{Simple value return*/
		/*!Return the Data<double> Energy*/
		Data<double> const& get_energy() const { return MCS_->get_energy(); }
		/*!Returns param_*/
		Vector<double> const& get_param() const { return param_; }
		/*!Return the coupling strength*/
		Vector<double> const& get_J() const { return MCS_->get_J(); }
		/*!Returns true if MCS_ can be run by MonteCarlo*/
		bool is_created() const { return (MCS_.get() && !MCS_->get_status()); }
		/*}*/

		/*!Analyse and output in ios*/
		std::string analyse(unsigned int const& level, IOSystem* ios);

		/*{Static methods*/
		static bool sort_by_E(MCSim const& a, MCSim const& b);
		static bool sort_by_coupling(MCSim const& a, MCSim const& b);
		static unsigned int sort_for_merge(MCSim const& list, MCSim const& new_elem);
		static unsigned int sort_by_param_for_merge(Vector<double> const& a, Vector<double> const& b);
		/*}*/

	private:
		Vector<double> param_;			//!< variational parameters used to create the wavefunction
		std::unique_ptr<MCSystem> MCS_; //!< pointer on a MCSystem that can be measured via MonteCarlo
};
#endif
