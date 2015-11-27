#ifndef DEF_VMCSYSTEMATIC
#define DEF_VMCSYSTEMATIC

#include "VMCMinimization.hpp"

class VMCSystematic : public VMCMinimization{
	public:
		VMCSystematic(VMCMinimization const& m, Vector<double> const& param, Matrix<int> const& sym, unsigned int const& p1, unsigned int const& p2);
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
		Vector<double> param_;
		Matrix<int> sym_;
		unsigned int p1_;
		unsigned int p2_;

		void apply_symmetry();
};
#endif
