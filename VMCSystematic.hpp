#ifndef DEF_VMCSYSTEMATIC
#define DEF_VMCSYSTEMATIC

#include "VMCMinimization.hpp"

class VMCSystematic : public VMCMinimization{
	public:
		VMCSystematic(VMCMinimization const& m);
		/*!Default destructor*/
		virtual ~VMCSystematic() = default;
		/*{Forbidden*/
		VMCSystematic() = delete;
		VMCSystematic(VMCSystematic const&) = delete;
		VMCSystematic(VMCSystematic&&) = delete;
		VMCSystematic& operator=(VMCSystematic) = delete;
		/*}*/

		void run(int const& nobs, double const& dE, unsigned int const& maxiter);
		void plot();
		void test();

	private:
		unsigned int maxiter_;
		int nobs_;
		double dE_;

		bool go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSystematic::*f)(Vector<double>*, Vector<unsigned int> const&));
		void evaluate(Vector<double>* x, Vector<unsigned int> const& idx);
};
#endif
