#ifndef DEF_VMCMINIMIZATION
#define DEF_VMCMINIMIZATION

#include "List.hpp"
#include "MCSim.hpp"

class VMCMinimization{
	public:
		VMCMinimization(Parseur& P);
		/*!Default destructor*/
		virtual ~VMCMinimization() = default;
		/*{Forbidden*/
		VMCMinimization() = delete;
		VMCMinimization(VMCMinimization const&) = delete;
		VMCMinimization(VMCMinimization&&) = delete;
		VMCMinimization& operator=(VMCMinimization) = delete;
		/*}*/

		void complete_analysis(double const& converged_criterion);
		void save() const;

	protected:
		List<MCSim> all_results_;
		Container system_param_;
		RST pso_info_;
		unsigned int Nfreedom_;
		unsigned int tmax_;
		double min_;
		double max_;
		double dx_;

		std::string get_filename() const { return basename_+"_"+time_; }
		void set_time() { time_ = Time().date(); }

	private:
		std::string basename_;
		std::string time_;
};
#endif

