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

		virtual void set_ps(unsigned int const& i, Vector<double> const& ps);
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

				std::string  wf_;
				unsigned int N_;
				unsigned int m_;
				unsigned int n_;
				int          bc_;
				unsigned int Nfreedom_;
				unsigned int tmax_;
				Vector<double>* ps_; //<! parameter space
				List<MCSim> all_results_;
				Container system_param_;
				RST pso_info_;
		};

	protected:
		IOFiles* out_;
		std::shared_ptr<Minimization> m_;

		void set_time(){ time_ = Time().date(); }
		std::string get_filename() const { return basename_+"_"+time_; }

		/*!Real call to GenericSystem+MonteCarlo*/
		std::shared_ptr<MCSim> evaluate(Vector<double> const& param);
};
#endif
