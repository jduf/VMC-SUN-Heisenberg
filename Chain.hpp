#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "CreateSystem.hpp"

class Chain: public CreateSystem<double>{
	public:
		Chain(Parseur& P);
		~Chain();

	private:
		/*!Compute the hopping matrix for a chain*/
		void compute_T();
		void save(std::string filename);
};

#endif
