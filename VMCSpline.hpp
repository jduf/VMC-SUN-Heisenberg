#ifndef DEF_VMCSPLINE
#define DEF_VMCSPLINE

#include "PSpline.hpp"
#include "VMCMinimization.hpp"

class VMCSpline : public VMCMinimization {
	public:
		VMCSpline(Parseur& P);
		/*!Default destructor*/
		virtual ~VMCSpline() = default;
		/*{Forbidden*/
		VMCSpline() = delete;
		VMCSpline(VMCSpline const&) = delete;
		VMCSpline(VMCSpline&&) = delete;
		VMCSpline& operator=(VMCSpline) = delete;
		/*}*/

		void run();
		void plot();

	private:
		PSpline pspline_;

		void compute_border();
		void set_param(Vector<double>& param, Vector<unsigned int>& idx);
};
#endif
