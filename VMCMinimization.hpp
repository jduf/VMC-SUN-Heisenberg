#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "Minimization.hpp"

class VMCMinimization{
	public:
		VMCMinimization(Minimization& m, std::string const& prefix);
		/*!Default destructor*/
		virtual ~VMCMinimization() = default;
		/*{Forbidden*/
		VMCMinimization() = delete;
		VMCMinimization(VMCMinimization const&) = delete;
		VMCMinimization(VMCMinimization&&) = delete;
		VMCMinimization& operator=(VMCMinimization const&) = delete;
		/*}*/

		virtual void set_ps(unsigned int const& i, Vector<double> const& ps)=0;
		void move(VMCMinimization* min);
		void refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax);
		void complete_analysis(double const& convergence_criterion);
		void save() const;

	protected:
		Minimization& m_;

		std::string get_filename() const { return basename_+"_"+time_; }
		void set_time() { time_ = Time().date(); }

		std::shared_ptr<MCSim> compute_vmc(Vector<double> const& param);

	private:
		std::string basename_;
		std::string time_;
};
#endif
