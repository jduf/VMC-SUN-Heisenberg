#ifndef DEF_VMCSPLINE
#define DEF_VMCSPLINE

#include "PSpline.hpp"
#include "VMCMinimization.hpp"

class VMCSpline : public VMCMinimization{
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

	private:
		PSpline pspline_;
		double min_;
		double max_;
		double dx_;

		void compute_border(double const& border, unsigned int const& dir);
		void set_param(Vector<double>& param);
};
#endif

