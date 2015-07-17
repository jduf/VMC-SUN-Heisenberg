#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "List.hpp"
#include "MCSim.hpp"

class VMCMinimization{
	public:
		VMCMinimization(Parseur& P);
		VMCMinimization(VMCMinimization const& m, std::string const& prefix);
		/*!Default destructor*/
		virtual ~VMCMinimization() = default;
		/*{Forbidden*/
		VMCMinimization() = delete;
		VMCMinimization(VMCMinimization const&) = delete;
		VMCMinimization(VMCMinimization&&) = delete;
		VMCMinimization& operator=(VMCMinimization const&) = delete;
		/*}*/

		void refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax);
		void complete_analysis(double const& convergence_criterion);
		void save() const;
		void save_best(unsigned int const& nsave);
		void plot() const;

		virtual void print() const;

	private:
		std::string time_;
		std::string basename_;

		class Minimization{
			public:
				Minimization(Parseur& P);
				/*!Default destructor*/
				virtual ~Minimization();
				/*{Forbidden*/
				Minimization() = delete;
				Minimization(Minimization const&) = delete;
				Minimization(Minimization&&) = delete;
				Minimization& operator=(Minimization const&) = delete;
				/*}*/
				bool within_limit(Vector<double> const& x);
				void save(IOFiles& out) const;

				List<MCSim> samples_list_;
				Container system_param_;
				Vector<double> J_; 
				RST pso_info_;
				double effective_time_ = 0.0;
				unsigned int tmax_     = 0;
				std::string  wf_       ="";
				unsigned int N_        = 0;
				unsigned int m_        = 0;
				unsigned int n_        = 0;
				int          bc_       = 0;
				unsigned int Nfreedom_ = 0;
				double       ps_size_  = 0;   //!< parameter space size
				Vector<double>* ps_    = NULL;//!< parameter space
		};

	protected:
		IOFiles* out_;
		std::shared_ptr<Minimization> m_;

		void set_time(){ time_ = Time().date("-"); }
		std::string get_filename() const { return time_+"_"+basename_; }

		/*!Real call to GenericSystem+MonteCarlo*/
		std::shared_ptr<MCSim> evaluate(Vector<double> const& param);
};
#endif
