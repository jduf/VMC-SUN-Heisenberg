#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "List.hpp"
#include "MCSim.hpp"

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
		void set_tmax(unsigned int const& tmax){ m_->tmax_ = tmax; };

		void refine();
		void refine(double const& E, double const& dE);
		void complete_analysis(double const& convergence_criterion);
		void save() const;
		void find_minima(unsigned int const& max_n_minima, List<MCSim>& list_min, Vector<double>& best_param, double& E_range) const;
		void find_and_run_minima(unsigned int const& max_n_minima);
		void find_save_and_plot_minima(unsigned int const& max_n_minima, IOFiles& w, std::string path="", std::string filename="") const;

		virtual void print() const;
		bool ready(){ return m_.get(); }
		RST& get_header(){ return m_->info_; }

	private:
		mutable std::string time_;
		std::string path_;
		std::string basename_;
		std::string prefix_;

		class Minimization{
			public:
				Minimization()=default;
				/*!Default destructor*/
				virtual ~Minimization();
				/*{Forbidden*/
				Minimization(Minimization const&) = delete;
				Minimization(Minimization&&) = delete;
				Minimization& operator=(Minimization const&) = delete;
				/*}*/
				void set(Parseur& P, std::string& path, std::string& basename);

				void create(Parseur& P, std::string& path, std::string& basename);
				std::string load(IOFiles& in, std::string& path, std::string& basename);
				void set_phase_space(Parseur const& P);

				bool within_limit(Vector<double> const& x);
				void save(IOFiles& out) const;

				List<MCSim> samples_list_;
				Container system_param_;
				RST info_;
				System* s_             = NULL;
				unsigned int dof_ 	   = 0;
				Vector<double>* ps_    = NULL;//!< parameter space
				double ps_size_        = 0.0; //!< parameter space size
				double effective_time_ = 0.0;
				unsigned int tmax_     = 0;
				Vector<double> J_;
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
		std::shared_ptr<MCSim> evaluate(Vector<double> const& param, unsigned int const& obs=0);
};
#endif
