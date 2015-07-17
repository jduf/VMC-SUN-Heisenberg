#ifndef DEF_VMCINTERPOLATION
#define DEF_VMCINTERPOLATION

#include "Interpolation.hpp"
#include "VMCMinimization.hpp"

class VMCInterpolation : public VMCMinimization{
	public:
		VMCInterpolation(VMCMinimization const& m);
		/*!Default destructor*/
		virtual ~VMCInterpolation() = default;
		/*{Forbidden*/
		VMCInterpolation() = delete;
		VMCInterpolation(VMCInterpolation const&) = delete;
		VMCInterpolation(VMCInterpolation&&) = delete;
		VMCInterpolation& operator=(VMCInterpolation) = delete;
		/*}*/

		void init();
		void run(unsigned int const& explore_around_minima);

		void plot();
		void print();

	private:
		Interpolation interp_;
		std::vector<Vector<unsigned int> > list_min_idx_;

		/*!Uses interpolation to find some/all minima*/
		void search_minima();
		/*!Iterative method that goes through all parameter space*/
		bool go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCInterpolation::*f)(Vector<double>*, Vector<unsigned int> const&));
		/*!Adds entry to list_min_idx_ if it's supposed to be a minima*/
		void select_if_min(Vector<double>* x, Vector<unsigned int> const& idx);
		/*!Saves (x,extrapolate(x)) in out_*/
		void save_interp_data(Vector<double>* x, Vector<unsigned int> const& idx);

		/*!Creates param via x and idx then call VMCMinimization::evaluate(param)*/
		void evaluate(Vector<double>* x, Vector<unsigned int> const& idx);
};
#endif
