#ifndef DEF_BINNING
#define DEF_BINNING

#include "Vector.hpp"
#include "Gnuplot.hpp"

class Binning{
	public:
		/*!Default constructor*/
		Binning(); 
		/*!Other constructor*/
		Binning(unsigned int B, unsigned int b); 
		/*!Default Simple destructor*/
		~Binning();
		/*Set the whole class*/
		void set();

		/*Add sample to the bins*/
		void add_sample(double const& x);
		/*True if the simulation is converged*/
		bool is_converged(double const& tol);

		/*Simple accessors to get the mean value*/
		double get_mean() const {return x_;}
		/*Simple accessors to get the variance*/
		double get_variance() const {return dx_;}

		/*Store the variance for each bin in a file*/
		void log(std::string filename){log_=new Write(filename);}
		/*Plot the variance in function of 1/l*/
		void plot();
		
		void check();

	private:
		double x_;		//!< mean
		double dx_;		//!< variance

		unsigned int const B_;	//!< Minimum number of biggest bins needed to compute variance
		unsigned int const b_;	//!< l_+b_ rank of the biggest bin (b !> 30)
		unsigned int l_;		//!< rank of the "smallest" bin 
		unsigned int DPL_;		//!< 2^l_ maximum number of element in each bin of rank l_
		unsigned int dpl_;		//!< current number of element in each bin of rank l_
		unsigned int log_iter_;

		Vector<unsigned int> Ml_;//!< Number bins of rank l : Ml = M0/2^l
		Vector<double> m_bin_;	//!< Mean of the Binnings
		Vector<double>*bin_;	//!< Binnings

		bool recompute_xdx_usefull_;

		Write* log_;

		/*!Recursive method that add samples in the different bins*/
		void add_bin(unsigned int l, double a, double b);
};
#endif
