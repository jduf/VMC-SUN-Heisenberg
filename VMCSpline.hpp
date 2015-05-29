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
		void run();
		void plot();

	private:
		PSpline pspline_;
		Vector<double>* border_;
		std::vector<Vector<unsigned int> > min_idx_;
		IOFiles* out_;

		bool go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSpline::*f)(Vector<double>*, Vector<unsigned int> const&));
		
		void run_if_min(Vector<double>* x, Vector<unsigned int> const& idx);
		void save_spline_data(Vector<double>* x, Vector<unsigned int> const& idx);

		//unsigned int VMCSpline::compute_border_size(Vector<unsigned int> const& edges_length){
		//unsigned int N(edges_length.size());
		//bool hypercube(true);
		//for(unsigned int i(1);i<N;i++){
		//if(edges_length(0) != edges_length(i)){ 
		//i = N;
		//hypercube = false;
		//}
		//}
		//unsigned int s(0);
		//if(hypercube){
		//for(unsigned int i(1);i<N+1;i++){
		//s += hypercube_element(N-i,N)*(i%2?1:-1)*pow(10,N-i);
		//}
		//} else {
		//for(unsigned int i(1);i<N+1;i++){
		//Vector<unsigned int> coeff(my::comb(N,N-i,edges_length));
		//if(i%2){ s += coeff.sum()*hypercube_element(N-i,N)/coeff.size(); }
		//else{ s -= coeff.sum()*hypercube_element(N-i,N)/coeff.size(); }
		//}
		//}
		//return s;
		//}

		//unsigned int VMCSpline::hypercube_element(int m, int n){
		//if(m==0 && n==0){ return 1; }
		//if(m<0 || n<0 || n<m){ return 0; }
		//return 2*hypercube_element(m,n-1)+hypercube_element(m-1,n-1);
		//}
};
#endif
