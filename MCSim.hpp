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

		/*!Sets S_ to a new MCSystem created via C*/
		void create_S(System const* const s);
		/*!Sets S_ to a copy obtained via MCSystem::clone() run on S*/
		void copy_S(std::unique_ptr<MCSystem> const& S);

		Vector<double> const& get_param() const { return param_; }
		std::unique_ptr<MCSystem> const& get_S() const { return S_; }

		void set_observable(bool all);
		void free_memory();

		bool is_created() const { return (S_.get() && !S_->get_status()); }
		void run(unsigned int const& thermalization_steps, unsigned int const& tmax);
		bool check_conv(double const& convergence_criterion);
		void complete_analysis(double const& convergence_criterion);

		void write(IOFiles& w) const;
		/*!Save the result in a single file (wavefunction parameters, observables)*/
		void save(IOFiles& w) const;
		void print() const { std::cout<<param_<<" "<<S_->get_energy(); }

		static bool sort_by_E(MCSim const& a, MCSim const& b);
		static unsigned int sort_by_param_for_merge(MCSim const& list, MCSim const& new_elem);
		static void merge(MCSim& list, MCSim& new_elem);

	private:
		Vector<double> param_;
		std::unique_ptr<MCSystem> S_;
};
#endif
