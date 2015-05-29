#ifndef DEF_VMCSPLINE
#define DEF_VMCSPLINE

#include "PSpline.hpp"
#include "VMCMinimization.hpp"

class VMCSpline : public VMCMinimization {
	public:
		VMCSpline(Parseur& P);
		/*!Default destructor*/
		virtual ~VMCSpline();
		/*{Forbidden*/
		VMCSpline() = delete;
		VMCSpline(VMCSpline const&) = delete;
		VMCSpline(VMCSpline&&) = delete;
		VMCSpline& operator=(VMCSpline) = delete;
		/*}*/

		void init(bool border);
		void run(unsigned int const& explore_around_minima);
		void plot();

	private:
		PSpline pspline_;
		Vector<double>* border_;
		std::vector<Vector<unsigned int> > all_min_idx_;
		IOFiles* out_;

		bool go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSpline::*f)(Vector<double>*, Vector<unsigned int> const&));
		
		void run_if_min(Vector<double>* x, Vector<unsigned int> const& idx);
		void save_spline_data(Vector<double>* x, Vector<unsigned int> const& idx);
		void call_compute_vmc(Vector<double>* x, Vector<unsigned int> const& idx);
};
#endif
