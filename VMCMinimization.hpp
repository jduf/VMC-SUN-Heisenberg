#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "List.hpp"
#include "MCSim.hpp"

class VMCMinimization{
	public:
		VMCMinimization(Parseur& P);
		/*!Default destructor*/
		virtual ~VMCMinimization()=0;
		/*{Forbidden*/
		VMCMinimization() = delete;
		VMCMinimization(VMCMinimization const&) = delete;
		VMCMinimization(VMCMinimization&&) = delete;
		VMCMinimization& operator=(VMCMinimization const&) = delete;
		/*}*/

		void set_x(unsigned int const& i, Vector<double> const& x);
		void move(VMCMinimization* min);
		void refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax);
		void complete_analysis(double const& converged_criterion);
		void save() const;

	protected:
		unsigned int Nfreedom_;
		unsigned int tmax_;
		Vector<double>* x_;
		List<MCSim> all_results_;
		Container system_param_;
		RST pso_info_;

		std::string get_filename() const { return basename_+"_"+time_; }
		void set_time() { time_ = Time().date(); }

		static bool sort_per_energy(MCSim const& a, MCSim const& b){ 
			return a.get_S()->get_energy().get_x()<b.get_S()->get_energy().get_x();
		};

	private:
		std::string basename_;
		std::string time_;
};
#endif

