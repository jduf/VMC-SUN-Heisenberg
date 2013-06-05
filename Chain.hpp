#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "CreateSystem.hpp"
#include "Parseur.hpp"

class Chain: public CreateSystem<double>{
	public:
		Chain(Parseur& P);
		~Chain();

	private:
		/*!Compute the hopping matrix for a chain*/
		void compute_T();
};

#endif
