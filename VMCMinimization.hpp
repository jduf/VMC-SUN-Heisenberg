#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "List.hpp"
#include "MCSim.hpp"
#include "omp.h"

class VMCMinimization{
	public:
		VMCMinimization(Parseur& P);
		VMCMinimization(IOFiles& in);
		/*!Default destructor*/
		virtual ~VMCMinimization() = default;
		/*{Forbidden*/
		VMCMinimization() = delete;
		VMCMinimization(VMCMinimization const&) = delete;
		VMCMinimization(VMCMinimization&&) = delete;
		VMCMinimization& operator=(VMCMinimization const&) = delete;
		/*}*/

		void set_phase_space(Parseur const& P){ m_->set_phase_space(P); }
		void set_tmax(unsigned int const& tmax){ m_->tmax_ = tmax; }

		void refine();
		void refine(double const& E, double const& dE);
		void refine(unsigned int const& nmin, int const& nobs, double const& dE, unsigned int const& maxiter);

		void complete_analysis(double const& convergence_criterion);
		void save() const;
		void save_parameters(unsigned int nbest) const;
		void run_parameters(Parseur& P);

		double find_minima(unsigned int const& max_local_minima, double const& range, List<MCSim>& sorted_list, List<MCSim>& list_min) const;
		void find_and_run_minima(unsigned int const& max_samples, int const& nobs, double const& dE);
		void find_save_and_plot_minima(unsigned int const& max_samples, IOFiles& w, std::string path="", std::string filename="") const;
		void explore_around_minima(unsigned int const& max_local_minima, int const& nobs, double const& dE, double const& dx);

		void improve_bad_samples(double const& dE);
		void clean();

		bool ready() const { return m_.get(); }
		RST& get_header(){ return m_->info_; }

	private:
		mutable std::string time_;
		std::string path_;
		std::string basename_;
		std::string prefix_;

		class Minimization{
			public:
				Minimization() = default;
				/*!Default destructor*/
				virtual ~Minimization();
				/*{Forbidden*/
				Minimization(Minimization const&) = delete;
				Minimization(Minimization&&) = delete;
				Minimization& operator=(Minimization const&) = delete;
				/*}*/
				void set(Parseur& P, std::string& path, std::string& basename);

				void create(Parseur& P, std::string& path, std::string& basename);
				void load(IOFiles& in, std::string& path, std::string& basename);
				void set_phase_space(Parseur const& P);

				bool within_limit(Vector<double> const& x) const;
				void save(IOFiles& out, bool const& all) const;

				List<MCSim> samples_;
				Container system_param_;
				RST info_;
				System* s_             = NULL;
				unsigned int dof_ 	   = 0;
				Vector<double>* ps_    = NULL;//!< parameter space
				double ps_size_        = 0.0; //!< parameter space size
				double effective_time_ = 0.0;
				unsigned int tmax_     = 0;
				std::vector<Observable> obs_;
		};

	protected:
		VMCMinimization(VMCMinimization const& m, std::string const& prefix);

		IOFiles* out_;
		std::shared_ptr<Minimization> m_;

		std::string const& get_path() const { return path_; }
		std::string get_filename() const { return time_+"_"+prefix_+basename_; }
		void set_time() const { time_ = Time().date("-"); }

		/*!Real call to the MonteCarlo evaluation via MCSim*/
		std::shared_ptr<MCSim> evaluate(Vector<double> const& param, int const& obs);
		void evaluate_until_precision(Vector<double> const& param, int const& nobs, double const& dE, unsigned int const& maxiter);
};
#endif
