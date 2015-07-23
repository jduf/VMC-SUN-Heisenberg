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

		void set_phase_space(Parseur& P){ m_->set_phase_space(P); }

		void refine();
		void refine(double const& E, double const& dE);
		void complete_analysis(double const& convergence_criterion);
		void save() const;
		void save_best(unsigned int const& nsave);
		void plot() const;

		virtual void print() const;
		bool ready(){ return m_.get(); }

	private:
		std::string time_;
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
				std::string set(Parseur& P);
				void set_phase_space(Parseur& P);

				bool within_limit(Vector<double> const& x);
				void save(IOFiles& out) const;

				List<MCSim> samples_list_;
				Container system_param_;
				RST pso_info_;
				System* s_             = NULL;
				unsigned int Nfreedom_ = 0;
				Vector<double>* ps_    = NULL;//!< parameter space
				double ps_size_        = 0;   //!< parameter space size
				double effective_time_ = 0.0;
				unsigned int tmax_     = 0;
		};

	protected:
		IOFiles* out_;
		std::shared_ptr<Minimization> m_;

		void set_time(){ time_ = Time().date("-"); }
		std::string get_filename() const { return time_+"_"+prefix_+basename_; }

		/*!Real call to GenericSystem+MonteCarlo*/
		std::shared_ptr<MCSim> evaluate(Vector<double> const& param);
};
#endif
