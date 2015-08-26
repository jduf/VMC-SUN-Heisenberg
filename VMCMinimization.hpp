#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "List.hpp"
#include "MCSim.hpp"
#include "Interpolation.hpp"

class VMCMinimization{
	public:
		VMCMinimization(Parseur& P);
		VMCMinimization(VMCMinimization const& m, std::string const& prefix);
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
		void set_time(){ time_ = Time().date("-"); }

		void refine();
		void refine(double const& E, double const& dE);
		void complete_analysis(double const& convergence_criterion);
		void save() const;
		void find_minima(unsigned int const& max_n_minima, List<MCSim>& list_min, Vector<double>& param, double& E_range, Interpolation<double>* interp_Er=NULL) const;
		void find_and_run_minima(unsigned int const& max_n_minima);
		void find_save_and_plot_minima(unsigned int const& max_n_minima, IOFiles& w, std::string path="", std::string filename="") const;

		virtual void print() const;
		bool ready(){ return m_.get(); }

	private:
		std::string time_;
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
				RST pso_info_;
				System* s_             = NULL;
				unsigned int dof_ 	   = 0;
				Vector<double>* ps_    = NULL;//!< parameter space
				double ps_size_        = 0;   //!< parameter space size
				double effective_time_ = 0.0;
				unsigned int tmax_     = 0;
		};

	protected:
		IOFiles* out_;
		std::shared_ptr<Minimization> m_;

		std::string const& get_path() const { return path_; }
		std::string get_filename() const { return time_+"_"+prefix_+basename_; }

		/*!Real call to the MonteCarlo evaluation via MCSim*/
		std::shared_ptr<MCSim> evaluate(Vector<double> const& param, unsigned int const& obs=0);
};
#endif
