#ifndef DEF_VMCSPLINE
#define DEF_VMCSPLINE

#include "PSpline.hpp"
#include "VMCMinimization.hpp"

class VMCSpline : public VMCMinimization{
	public:
		VMCSpline(VMCMinimization const& m);
		/*!Default destructor*/
		virtual ~VMCSpline() = default;
		/*{Forbidden*/
		VMCSpline() = delete;
		VMCSpline(VMCSpline const&) = delete;
		VMCSpline(VMCSpline&&) = delete;
		VMCSpline& operator=(VMCSpline) = delete;
		/*}*/

		void init();
		void run(unsigned int const& explore_around_minima);

		void plot();
		void print();

	private:
		PSpline pspline_;
		std::vector<Vector<unsigned int> > list_min_idx_;

		/*!Uses a PSpline to find all minima*/
		void search_minima();
		/*!Iterative method that goes through all parameter space*/
		bool go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSpline::*f)(Vector<double>*, Vector<unsigned int> const&));
		/*!Adds entry to list_min_idx_ if it's supposed to be a minima*/
		void select_if_min(Vector<double>* x, Vector<unsigned int> const& idx);
		/*!Saves (x,extrapolate(x)) in out_*/
		void save_spline_data(Vector<double>* x, Vector<unsigned int> const& idx);

		/*!Creates param via x and idx then call VMCMinimization::evaluate(param)*/
		void evaluate(Vector<double>* x, Vector<unsigned int> const& idx);
};
#endif
