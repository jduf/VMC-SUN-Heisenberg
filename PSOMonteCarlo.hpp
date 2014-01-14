#ifndef DEF_PSOMONTECARLO
#define DEF_PSOMONTECARLO

#include "MonteCarlo.hpp"
#include "SquareJastrow.hpp"
#include "TriangleJastrow.hpp"
#include "PSO.hpp"

class PSOMonteCarlo : public PSO {
	public:
		PSOMonteCarlo(Parseur& P,unsigned int Nfreedom);
		virtual ~PSOMonteCarlo();

		void save();

	private:
		virtual double run(Vector<double> const& x);
		void copy_input(Container const& input, Container& new_input);
		void get_system_properties(Parseur& P, Container& c);
		Matrix<double> create_nu(Vector<double> const& x);

		Container input_;
		Container param_;
		Write results_;
		std::string wf_;
		unsigned int z_;
};
#endif

