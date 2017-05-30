#ifndef DEF_VMCSYSTEMATIC
#define DEF_VMCSYSTEMATIC

#include "VMCMinimization.hpp"

class VMCSystematic : public VMCMinimization{
	public:
		VMCSystematic(VMCMinimization const& m);
		VMCSystematic(IOFiles& in);
		/*!Default destructor*/
		virtual ~VMCSystematic() = default;
		/*{Forbidden*/
		VMCSystematic() = delete;
		VMCSystematic(VMCSystematic const&) = delete;
		VMCSystematic(VMCSystematic&&) = delete;
		VMCSystematic& operator=(VMCSystematic) = delete;
		/*}*/

		void run(double const& dEoE, unsigned int const& tmax, unsigned int const& maxiter, std::string const& save_in);
		void plot();
		void analyse(std::string const& path, std::string const& filename, List<MCSim>& keep) const;

	private:
		unsigned int maxiter_= 0;
		unsigned int tmax_	 = 0;
		double dEoE_ 		 = 0.0;

		bool go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSystematic::*f)(Vector<double>*, Vector<unsigned int> const&));
		void evaluate(Vector<double>* x, Vector<unsigned int> const& idx);
};
#endif
