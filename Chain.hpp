#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "CreateSystem.hpp"

class Chain: public CreateSystem<double>{
	public:
		Chain(Parseur& P);
		~Chain();

	private:
		void compute_H();
		void compute_EVec();
		void save(std::string filename);
};

#endif
