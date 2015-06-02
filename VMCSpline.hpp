#ifndef DEF_VMCSPLINE
#define DEF_VMCSPLINE

#include "PSpline.hpp"
#include "VMCMinimization.hpp"

class VMCSpline : public VMCMinimization {
	public:
		VMCSpline(Minimization& m);
		/*!Default destructor*/
		virtual ~VMCSpline() = default;
		/*{Forbidden*/
		VMCSpline() = delete;
		VMCSpline(VMCSpline const&) = delete;
		VMCSpline(VMCSpline&&) = delete;
		VMCSpline& operator=(VMCSpline) = delete;
		/*}*/

		void init();
		void set_ps(unsigned int const& i, Vector<double> const& ps);
		void run(unsigned int const& explore_around_minima);
		void plot();
		void print() const;

	private:
		PSpline pspline_;
		std::vector<Vector<unsigned int> > all_min_idx_;
		IOFiles* out_;

		void search_minima();
		bool go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSpline::*f)(Vector<double>*, Vector<unsigned int> const&));
		
		void select_if_min(Vector<double>* x, Vector<unsigned int> const& idx);
		void save_spline_data(Vector<double>* x, Vector<unsigned int> const& idx);
		void call_compute_vmc(Vector<double>* x, Vector<unsigned int> const& idx);
};
#endif
