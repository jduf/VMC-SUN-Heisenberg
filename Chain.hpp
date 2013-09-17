#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "CreateSystem.hpp"
#include "Gnuplot.hpp"

class Chain: public CreateSystem<double>{
	public:
		Chain(Parseur& P);
		~Chain();

	private:
		void compute_H();
		void compute_P();
		void compute_T();
		void compute_EVec();
		void compute_spectrum();
		void save(std::string filename);

		Matrix<double> Px;
};
#endif
