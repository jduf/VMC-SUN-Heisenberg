#ifndef DEF_VMCSYSTEMATIC
#define DEF_VMCSYSTEMATIC

#include "VMCMinimization.hpp"

class VMCSystematic : public VMCMinimization{
	public:
		VMCSystematic(VMCMinimization const& m, Parseur& P);
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
};
#endif
