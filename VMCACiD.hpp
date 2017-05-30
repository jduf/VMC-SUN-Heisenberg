#ifndef ROSENBOCKACID
#define ROSENBOCKACID

#include "ACiD.hpp"
#include "VMCMinimization.hpp"

class VMCACiD: public VMCMinimization,  public ACiD{
	public:
		/*!Constructor*/
		VMCACiD(VMCMinimization const& min, Vector<unsigned int> const& d, Vector<double> const& param);
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
		void run(double const& dEoE, unsigned int const& tmax, unsigned int const& maxiter, unsigned int const& maxsteps, std::string const& save_in);
		bool keepon(bool const& improve_overall) const;

		void display_param_and_xmean(Vector<double> const& param) const;
		void test();
		void save(std::string dirname) const;

	private:
		Matrix<int> idx_;
		unsigned int maxiter_ = 0;
		unsigned int tmax_ = 0;
		double dEoE_;
};
#endif
