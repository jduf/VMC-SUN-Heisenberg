#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "MCSim.hpp"
#include "omp.h"

class VMCMinimization{
	public:
		/*!Constructor that start a new minimisation from arguments given via Parseur*/
		VMCMinimization(Parseur& P);
		/*!Constructor that reads data from in and loads all m_->samples if loadall=true*/
		VMCMinimization(IOFiles& in, bool const& loadall, std::string const& prefix);
		/*!Default destructor*/
		virtual ~VMCMinimization() = default;
		/*{Forbidden*/
		VMCMinimization() = delete;
		VMCMinimization(VMCMinimization const&) = delete;
		VMCMinimization(VMCMinimization&&) = delete;
		VMCMinimization& operator=(VMCMinimization const&) = delete;
		/*}*/

		void load(IOFiles& in, bool const& loadall){ m_->load(in,path_,basename_,loadall); }
		void set_phase_space(Parseur const& P){ m_->set_phase_space(P); }
		void swap_phase_space(Vector<double>*& ps){ m_->swap_phase_space(ps); }
		void set_tmax(unsigned int const& tmax){ m_->tmax_ = tmax; }

		/*!Different ways to improve the measures of the best samples, saves the result*/
		void refine(std::string const& save_in);
		void refine(double const& E, double const& dEoE, std::string const& save_in);
		void refine(unsigned int const& nmin, Vector<unsigned int> const& which_obs, double const& dEoE, unsigned int const& maxiter, std::string const& save_in);
		void refine(unsigned int const& ttotal, Vector<unsigned int> const& which_obs, double const& dEoE, std::string const& save_in);

		void improve_bad_samples(double const& dEoE, std::string const& save_in);

		void complete_analysis(double const& convergence_criterion);
		virtual void save(std::string save_in) const;

		double find_minima(unsigned int const& max_local_minima, double const& range, List<MCSim>& sorted_list, List<MCSim>& list_min) const;
		void find_and_run_minima(unsigned int const& max_samples, Vector<unsigned int> const& which_obs, double const& dEoE, std::string const& save_in);
		void find_save_and_plot_minima(unsigned int const& max_samples, IOFiles& w, std::string path="", std::string filename="") const;
		void explore_around_minima(unsigned int const& max_local_minima, Vector<unsigned int> const& which_obs, double const& dEoE, double const& dx, std::string const& save_in);

		void save_parameters(unsigned int nbest, std::string save_in) const;
		void run_parameters(std::string const& parameters_file, Vector<unsigned int> const& which_obs, double const& dEoE, unsigned int const& maxiter, std::string const& save_in);

		RST& get_header(){ return m_->info_; }
		bool ready() const { return m_.get(); }
		std::string const& get_path() const { return path_; }
		std::string get_filename() const { return time_+"_"+prefix_+basename_; }
		unsigned int get_number_of_samples() const { return m_.get()?m_->samples_.size():0; }
		List<MCSim>::Node* get_best_target() const;

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
				void load(IOFiles& in, std::string& path, std::string& filename, bool const& loadall);
				bool set_phase_space(Parseur const& P);
				bool swap_phase_space(Vector<double>*& ps);

				bool within_limit(Vector<double> const& x) const;
				void save(IOFiles& out, bool const& all) const;

				List<MCSim> samples_;
				RST info_;
				System* s_             = NULL;
				unsigned int dof_ 	   = 0;   //!< number of variational parameters (degrees of freedom)
				Vector<double>* ps_    = NULL;//!< parameter space
				double ps_size_        = 0.0; //!< parameter space size
				double effective_time_ = 0.0;
				unsigned int tmax_     = 0;
				std::vector<Observable> obs_;
		};

	protected:
		VMCMinimization(VMCMinimization const& m, std::string const& prefix);
		void load_filenames(IOFiles& in){ in>>path_>>basename_; }

		IOFiles* out_;
		std::shared_ptr<Minimization> m_;
		unsigned int progress_;
		unsigned int total_eval_;

		void set_time() const { time_ = Time().date("-"); }
		std::string const& get_time() const { return time_; }

		/*!Real call to the MonteCarlo evaluation via MCSim*/
		std::shared_ptr<MCSim> evaluate(Vector<double> const& param, Vector<unsigned int> const& which_obs);
		std::shared_ptr<MCSim> evaluate_until_precision(Vector<double> const& param, Vector<unsigned int> const& which_obs, double const& dEoE, unsigned int const& maxiter);
		void save(IOFiles& out) const;
};
#endif
