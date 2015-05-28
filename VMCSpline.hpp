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

		void init(bool border);
		void run();
		void plot();

	private:
		PSpline pspline_;
		std::vector<Vector<double> > border_;
		std::vector<Vector<double> > parameter_space_;
		std::vector<Vector<unsigned int> > parameter_idx_;

		bool set_parameter_space(Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, bool const& border);
};
#endif
