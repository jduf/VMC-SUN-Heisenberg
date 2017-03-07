#ifndef ROSENBOCKACID
#define ROSENBOCKACID

#include "ACiD.hpp"
#include "VMCMinimization.hpp"

class VMCACiD: public VMCMinimization,  public ACiD{
	public:
		/*!Constructor*/
		VMCACiD(VMCMinimization const& min, Vector<unsigned int> const& d);
		/*!Constructor that reads from file*/
		VMCACiD(VMCMinimization const& min, IOFiles& in_ACiD);
		/*!Default destructor*/
		~VMCACiD()=default;
		/*{Forbidden*/
		VMCACiD() = delete;
		VMCACiD(VMCACiD const&) = delete;
		VMCACiD(VMCACiD&&) = delete;
		VMCACiD& operator=(VMCACiD);
		/*}*/

		double function(Vector<double> const& x);
		void run(unsigned int const& maxsteps, unsigned int const& maxiter, double const& dEoE);
		void save(std::string dirname);

	private:
		Matrix<int> idx_;
		double dEoE_;
		unsigned int maxiter_;
};
#endif
