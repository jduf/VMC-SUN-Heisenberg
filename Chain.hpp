#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "CreateSystem.hpp"

class Chain: public CreateSystem<double>{
	public:
		Chain(Parseur& P);
		~Chain();

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();
		void save(std::string filename);

		Vector<unsigned int> get_neighbourg(unsigned int i);
};
#endif
